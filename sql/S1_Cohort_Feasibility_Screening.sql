-- ============================================================
-- MIMIC-IV Feasibility Screening: Cancer + KDIGO >= 2 AKI + RRT Initiation Timing (Counts only)
-- Syntax: Google BigQuery (Official PhysioNet BigQuery mirror)
-- Before running, please ensure:
--   1) Replace the dataset names below with the version you have access to, e.g.,
--      physionet-data.mimiciv_3_1_hosp / _icu / _derived
--      (If the 'derived' dataset is unavailable, generate kdigo_stages, rrt, icustay_detail using the mimic-code repository first)
--   2) Perform a random sample check to verify that 'dialysis_active' in derived.rrt matches 'chartevents'
-- ============================================================

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

-- Number of creatinine measurements per patient within 7 days after t0 (and during ICU stay)
, creat AS (
  SELECT f.stay_id, COUNT(*) AS n_cr
  FROM final f
  JOIN icu i USING (stay_id)
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON l.hadm_id = f.hadm_id
   AND l.itemid = 50912                    -- Serum creatinine
   AND l.valuenum IS NOT NULL
   AND l.charttime BETWEEN f.t0
       AND LEAST(DATETIME_ADD(f.t0, INTERVAL 7 DAY), i.icu_outtime)
  GROUP BY f.stay_id
),
joined AS (
  SELECT f.stay_id, f.solid, f.heme, f.hospital_expire_flag,
         COALESCE(c.n_cr, 0) AS n_cr
  FROM final f
  LEFT JOIN creat c USING (stay_id)
)

-- ============================================================
-- Output: Whether the number of creatinine measurements is sufficient for trajectory clustering
-- ============================================================
SELECT
  COUNT(*)                                   AS total,
  COUNTIF(n_cr >= 3)                         AS n_ge3,
  COUNTIF(n_cr >= 5)                         AS n_ge5,
  COUNTIF(n_cr >= 7)                         AS n_ge7,
  APPROX_QUANTILES(n_cr, 2)[OFFSET(1)]       AS median_n,
  COUNTIF(n_cr >= 3 AND heme = 1)            AS ge3_any_heme,
  COUNTIF(n_cr >= 3 AND solid = 1 AND heme = 0) AS ge3_solid_only,
  COUNTIF(n_cr >= 3 AND hospital_expire_flag = 1) AS ge3_hosp_deaths
FROM joined;
