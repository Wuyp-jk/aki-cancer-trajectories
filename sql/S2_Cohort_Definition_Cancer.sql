-- ============================================================
-- 04a  Create Cohort Table (Run only once)
-- Purpose: Save your validated cohort (4,213 patients) as a table in your own project,
--          so that subsequent scripts (04b~04e) can query from this table without repeating the long SQL.
-- Notes:
--   1) The project name 'YOUR_PROJECT' is from your previous screenshot; replace it globally if different.
--   2) In the free sandbox mode, tables expire automatically after 60 days. Please save the exported CSVs locally.
--   3) The dataset location must be US (to match MIMIC).
-- ============================================================

CREATE SCHEMA IF NOT EXISTS `YOUR_PROJECT.aki_cancer`
OPTIONS (location = 'US');

CREATE OR REPLACE TABLE `YOUR_PROJECT.aki_cancer.cohort` AS
WITH
-- 1. Adults, first hospital admission, first ICU stay
icu AS (
  SELECT stay_id, subject_id, hadm_id, icu_intime, icu_outtime,
         admission_age, hospital_expire_flag
  FROM `physionet-data.mimiciv_3_1_derived.icustay_detail`
  WHERE admission_age >= 18
    AND first_hosp_stay
    AND first_icu_stay
),

-- 2. Cancer (ICD codes, hospital diagnosis level; Note: cannot distinguish "active" cancer, must be addressed in limitations)
cancer AS (
  SELECT hadm_id,
    MAX(CASE
          WHEN icd_version = 10 AND SUBSTR(icd_code,1,1) = 'C'
               AND SUBSTR(icd_code,1,3) BETWEEN 'C81' AND 'C96' THEN 1
          WHEN icd_version = 9
               AND SAFE_CAST(SUBSTR(icd_code,1,3) AS INT64) BETWEEN 200 AND 208 THEN 1
          ELSE 0 END) AS heme,
    MAX(CASE
          WHEN icd_version = 10 AND SUBSTR(icd_code,1,1) = 'C'
               AND SUBSTR(icd_code,1,3) BETWEEN 'C00' AND 'C80'
               AND SUBSTR(icd_code,1,3) != 'C44' THEN 1
          WHEN icd_version = 9
               AND SAFE_CAST(SUBSTR(icd_code,1,3) AS INT64) BETWEEN 140 AND 199
               AND SUBSTR(icd_code,1,3) != '173' THEN 1
          ELSE 0 END) AS solid
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  GROUP BY hadm_id
  HAVING heme = 1 OR solid = 1
),

-- 3. Pre-existing ESRD / Chronic dialysis (Exclusion)
esrd AS (
  SELECT DISTINCT hadm_id
  FROM `physionet-data.mimiciv_3_1_hosp.diagnoses_icd`
  WHERE (icd_version = 10 AND icd_code IN ('N186','Z992'))
     OR (icd_version = 9  AND icd_code IN ('5856','V451','V4511','V4512'))
),

-- 4. Time zero (t0): First occurrence of KDIGO >= 2 (Using only creatinine/urine output criteria, excluding CRRT in aki_stage, 
--    otherwise patients receiving RRT would be artificially classified as stage 3, causing a circular definition)
t0 AS (
  SELECT stay_id, MIN(charttime) AS t0
  FROM `physionet-data.mimiciv_3_1_derived.kdigo_stages`
  WHERE GREATEST(COALESCE(aki_stage_creat,0), COALESCE(aki_stage_uo,0)) >= 2
  GROUP BY stay_id
),

-- 5. Time of first RRT
rrt1 AS (
  SELECT stay_id, MIN(charttime) AS rrt_time
  FROM `physionet-data.mimiciv_3_1_derived.rrt`
  WHERE dialysis_active = 1
  GROUP BY stay_id
),

-- 6. Base cohort
cohort AS (
  SELECT i.stay_id, i.hadm_id, i.hospital_expire_flag,
         c.solid, c.heme, t.t0, r.rrt_time,
         DATETIME_DIFF(r.rrt_time, t.t0, HOUR) AS rrt_h
  FROM icu i
  JOIN cancer c USING (hadm_id)
  JOIN t0 t USING (stay_id)
  LEFT JOIN rrt1 r USING (stay_id)
  WHERE i.hadm_id NOT IN (SELECT hadm_id FROM esrd)
),

-- 7. Urgent indications: K > 6.0 or pH < 7.15 within 6 hours before t0 (Window can be adjusted according to your protocol)
urgent AS (
  SELECT DISTINCT co.stay_id
  FROM cohort co
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON l.hadm_id = co.hadm_id
   AND l.charttime BETWEEN DATETIME_SUB(co.t0, INTERVAL 6 HOUR) AND co.t0
  WHERE (l.itemid IN (50971, 52610) AND l.valuenum > 6.0)   -- Potassium
     OR (l.itemid = 50820 AND l.valuenum < 7.15)            -- Blood gas pH
),

-- 8. Final analytical cohort (Excluding: RRT before t0, presence of urgent indications)
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
  WHERE (co.rrt_time IS NULL OR co.rrt_time >= co.t0)       -- Exclude RRT before t0
    AND co.stay_id NOT IN (SELECT stay_id FROM urgent)
)

-- The following section creates a table with the cohort and baseline information (the cohort definition itself remains unchanged)
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

-- Validation: 'total' should be 4213, consistent with your previous run
SELECT COUNT(*) AS total,
       COUNT(DISTINCT stay_id) AS distinct_stays,
       COUNTIF(heme = 1) AS any_heme,
       COUNTIF(solid = 1 AND heme = 0) AS solid_only
FROM `YOUR_PROJECT.aki_cancer.cohort`;
