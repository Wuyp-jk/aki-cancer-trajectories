-- ============================================================
-- S3  Create "Non-cancer Control" Cohort Table (Run only once)
-- The only difference from 04a is that the cancer condition is inverted:
--   04a: Must have cancer diagnostic codes (JOIN cancer)
--   04f: Must NOT have any malignant tumor-related codes (NOT IN anycancer)
-- All other criteria (adults, first hospital admission, first ICU stay, t0 definition, 
-- ESRD exclusion, pre-t0 RRT exclusion, urgent indication exclusion) are identical to 04a.
-- Codes excluded from the control group: ICD-10 starting with 'C' (excluding C44) and Z85 (personal history of malignancy);
--                                       ICD-9 140-209 (excluding 173) and V10 (personal history of malignancy).
-- Therefore, patients with only C44 (non-melanoma skin cancer) are in neither group, which is intentional.
-- 'solid' / 'heme' are fixed at 0 in the control group to maintain the column structure for 04b and 04e.
-- ============================================================

CREATE OR REPLACE TABLE `YOUR_PROJECT.aki_cancer.cohort_ctrl` AS
WITH
icu AS (
  SELECT stay_id, subject_id, hadm_id, icu_intime, icu_outtime,
         admission_age, hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_derived.icustay_detail`
  WHERE admission_age >= 18
    AND first_hosp_stay
    AND first_icu_stay
),

anycancer AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 10 AND (SUBSTR(icd_code,1,1) = 'C' OR SUBSTR(icd_code,1,3) = 'Z85'))
     OR (icd_version = 9
         AND (SAFE_CAST(SUBSTR(icd_code,1,3) AS INT64) BETWEEN 140 AND 209
              OR SUBSTR(icd_code,1,3) = 'V10'))
),

esrd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 10 AND icd_code IN ('N186','Z992'))
     OR (icd_version = 9  AND icd_code IN ('5856','V451','V4511','V4512'))
),

t0 AS (
  SELECT stay_id, MIN(charttime) AS t0
  FROM `physionet-data.mimiciv_3_1_derived.kdigo_stages`
  WHERE GREATEST(COALESCE(aki_stage_creat,0), COALESCE(aki_stage_uo,0)) >= 2
  GROUP BY stay_id
),

rrt1 AS (
  SELECT stay_id, MIN(charttime) AS rrt_time
  FROM `physionet-data.mimiciv_3_1_derived.rrt`
  WHERE dialysis_active = 1
  GROUP BY stay_id
),

cohort AS (
  SELECT i.stay_id, i.hadm_id, i.hospital_expire_flag,
         0 AS solid, 0 AS heme, t.t0, r.rrt_time,
         DATETIME_DIFF(r.rrt_time, t.t0, HOUR) AS rrt_h
  FROM icu i
  JOIN t0 t USING (stay_id)
  LEFT JOIN rrt1 r USING (stay_id)
  WHERE i.hadm_id NOT IN (SELECT hadm_id FROM anycancer)
    AND i.hadm_id NOT IN (SELECT hadm_id FROM esrd)
),

urgent AS (
  SELECT DISTINCT co.stay_id
  FROM cohort co
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON l.hadm_id = co.hadm_id
   AND l.charttime BETWEEN DATETIME_SUB(co.t0, INTERVAL 6 HOUR) AND co.t0
  WHERE (l.itemid IN (50971, 52610) AND l.valuenum > 6.0)
     OR (l.itemid = 50820 AND l.valuenum < 7.15)
),

final AS (
  SELECT co.*,
    CASE
      WHEN co.rrt_h IS NULL                THEN 'never_rrt'
      WHEN co.rrt_h <= 12                  THEN 'rrt_within_12h'
      WHEN co.rrt_h <= 24                  THEN 'rrt_12_24h'
      WHEN co.rrt_h <= 72                  THEN 'rrt_24_72h'
      ELSE                                      'rrt_after_72h'
    END AS rrt_group
  FROM cohort co
  WHERE (co.rrt_time IS NULL OR co.rrt_time >= co.t0)
    AND co.stay_id NOT IN (SELECT stay_id FROM urgent)
)

SELECT
  f.stay_id, f.hadm_id, d.subject_id,
  d.gender, d.race, d.admission_age,
  f.solid, f.heme,
  f.t0,
  d.icu_intime, d.icu_outtime, d.admittime, d.dischtime, d.dod,
  d.los_icu, d.los_hospital,
  f.hospital_expire_flag,
  f.rrt_time, f.rrt_h, f.rrt_group
FROM final f
JOIN `physionet-data.mimiciv_3_1_derived.icustay_detail` d
  ON d.stay_id = f.stay_id;

-- Validation: 'total' is the control group size; 'overlap_with_cancer' must be 0
SELECT COUNT(*) AS total,
       COUNT(DISTINCT stay_id) AS distinct_stays,
       COUNTIF(stay_id IN (SELECT stay_id FROM `YOUR_PROJECT.aki_cancer.cohort`)) AS overlap_with_cancer
FROM `YOUR_PROJECT.aki_cancer.cohort_ctrl`;
