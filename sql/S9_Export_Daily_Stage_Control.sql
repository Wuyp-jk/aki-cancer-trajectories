-- ============================================================
-- 04g  Daily Maximum KDIGO Stage (Non-cancer Control Cohort)
-- Purpose: Export the highest KDIGO stage per day per patient (max 7 rows per patient).
-- Export file name: daily_stage_ctrl.csv
-- Note: Uses cohort_ctrl (non-cancer).
-- ============================================================
SELECT
  k.stay_id,
  LEAST(CAST(FLOOR(DATETIME_DIFF(k.charttime, c.t0, MINUTE) / 60.0 / 24) AS INT64), 6) AS day_n,
  MAX(CASE
        WHEN k.aki_stage_creat IS NULL AND k.aki_stage_uo IS NULL THEN NULL
        ELSE GREATEST(COALESCE(k.aki_stage_creat, 0), COALESCE(k.aki_stage_uo, 0))
      END) AS max_stage
FROM `YOUR_PROJECT.aki_cancer.cohort_ctrl` c
JOIN `physionet-data.mimiciv_3_1_derived.kdigo_stages` k
  ON k.stay_id = c.stay_id
WHERE DATETIME_DIFF(k.charttime, c.t0, MINUTE) BETWEEN 0 AND 168 * 60
GROUP BY k.stay_id, day_n
ORDER BY k.stay_id, day_n;
