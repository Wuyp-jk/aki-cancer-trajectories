-- ============================================================
-- 04e  Outcomes Table (One row per patient)
-- Export file name: outcomes.csv
-- Notes:
--   'dod' (date of death) is from MIMIC's post-discharge mortality records (up to ~1 year post-discharge).
--   'death_28d' / 'death_90d' are calculated from t0. If 'dod' is missing, the patient is considered alive during the follow-up period.
--   This assumption must be stated in the limitations section of the manuscript.
--   'last_creat' is the last creatinine measured before hospital discharge, used to assess renal recovery.
--   The specific definition of "recovery" will be determined by comparing it with 'creat_baseline' during analysis.
-- ============================================================

WITH last_cr AS (
  SELECT l.hadm_id,
         l.valuenum  AS last_creat,
         l.charttime AS last_creat_time
  FROM `physionet-data.mimiciv_3_1_hosp.labevents` l
  JOIN `YOUR_PROJECT.aki_cancer.cohort` c ON c.hadm_id = l.hadm_id
  WHERE l.itemid = 50912
    AND l.valuenum IS NOT NULL
    AND l.charttime >= c.t0
    AND l.charttime <= c.dischtime
  QUALIFY ROW_NUMBER() OVER (PARTITION BY l.hadm_id ORDER BY l.charttime DESC) = 1
)
SELECT
  c.stay_id, c.hadm_id,
  c.admittime, c.dischtime, c.t0, c.icu_outtime,
  a.deathtime,
  c.dod,
  c.hospital_expire_flag,
  c.los_icu, c.los_hospital,
  (c.rrt_time IS NOT NULL)                                   AS any_rrt,
  c.rrt_time, c.rrt_h, c.rrt_group,
  lc.last_creat, lc.last_creat_time,
  (c.dod IS NOT NULL AND DATE_DIFF(c.dod, DATE(c.t0), DAY) <= 28) AS death_28d,
  (c.dod IS NOT NULL AND DATE_DIFF(c.dod, DATE(c.t0), DAY) <= 90) AS death_90d
FROM `YOUR_PROJECT.aki_cancer.cohort` c
LEFT JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.hadm_id = c.hadm_id
LEFT JOIN last_cr lc ON lc.hadm_id = c.hadm_id
ORDER BY c.stay_id;
