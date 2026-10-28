-- ============================================================
-- 04d  KDIGO Time Series (Creatinine + Urine Output Rate, much denser than creatinine alone)
-- Purpose: Compensate for the issue of sparse creatinine measurements (median of only 3)
-- Time window: 24 hours before t0 to 7 days after t0 (capped at ICU discharge time)
-- Column description:
--   This query extracts all existing columns from the 'kdigo_stages' table for this version.
--   (Typically includes creat, urine output rates uo_rt_6hr/12hr/24hr, aki_stage_creat, aki_stage_uo, etc.)
--   This prevents errors caused by misspelled column names. Once exported, check the actual columns in R.
-- Export file name: kdigo_series.csv
-- ============================================================

SELECT
  c.stay_id,
  k.charttime,
  ROUND(DATETIME_DIFF(k.charttime, c.t0, MINUTE) / 60.0, 2) AS hours_from_t0,
  k.* EXCEPT (subject_id, hadm_id, stay_id, charttime)
FROM `YOUR_PROJECT.aki_cancer.cohort` c
JOIN `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
  ON k.stay_id = c.stay_id
 AND k.charttime BETWEEN DATETIME_SUB(c.t0, INTERVAL 24 HOUR)
     AND LEAST(DATETIME_ADD(c.t0, INTERVAL 7 DAY), c.icu_outtime)
ORDER BY c.stay_id, k.charttime;
