-- ============================================================
-- 04c  Raw Creatinine Sequence (One row per measurement)
-- Time window: 48 hours before t0 to 7 days after t0 (capped at ICU discharge time)
-- Export file name: creatinine_raw.csv
-- ============================================================

SELECT
  c.stay_id,
  l.charttime,
  ROUND(DATETIME_DIFF(l.charttime, c.t0, MINUTE) / 60.0, 2) AS hours_from_t0,
  l.valuenum AS creat
FROM `YOUR_PROJECT.aki_cancer.cohort` c
JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
  ON l.hadm_id = c.hadm_id
 AND l.itemid = 50912                       -- Serum creatinine
 AND l.valuenum IS NOT NULL
 AND l.charttime BETWEEN DATETIME_SUB(c.t0, INTERVAL 48 HOUR)
     AND LEAST(DATETIME_ADD(c.t0, INTERVAL 7 DAY), c.icu_outtime)
ORDER BY c.stay_id, l.charttime;
