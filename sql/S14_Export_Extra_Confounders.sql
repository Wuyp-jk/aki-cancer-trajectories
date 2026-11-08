-- ============================================================
-- S14  Extra baseline confounders for the residual-confounding sensitivity analysis
-- Purpose: admission type, first ICU unit, admitting service, calendar period and
--          code status at or before t0, for every patient of both cohort tables.
-- Needs:   the cohort tables created by S2 (cancer) and S3 (control).
-- Export:  save the result as extra_confounders.csv (one row per stay_id).
-- Notes:
--   * Code status uses the two MIMIC-IV chartevents items for "Code Status" (itemid 128 and 223758).
--     Only entries charted at or before t0 are used, so that orders written after AKI onset
--     (which lie on the path from cancer to outcome) do not enter the baseline covariates.
--     Patients with no entry before t0 get NULL (coded "not documented" in R).
--   * anchor_year_group is the only calendar information MIMIC-IV gives (3-year bands).
-- ============================================================
WITH coh AS (
  SELECT stay_id, subject_id, hadm_id, t0, 1 AS cancer
  FROM `YOUR_PROJECT.aki_cancer.cohort`
  UNION ALL
  SELECT stay_id, subject_id, hadm_id, t0, 0 AS cancer
  FROM `YOUR_PROJECT.aki_cancer.cohort_ctrl`
),
cs AS (
  SELECT c.stay_id,
         ce.value AS code_status_raw,
         ROW_NUMBER() OVER (PARTITION BY c.stay_id ORDER BY ce.charttime DESC) AS rn
  FROM coh c
  JOIN `physionet-data.mimiciv_3_1_icu.chartevents` ce
    ON ce.stay_id = c.stay_id
   AND ce.itemid IN (128, 223758)
   AND ce.charttime <= c.t0
),
svc AS (
  SELECT hadm_id, curr_service,
         ROW_NUMBER() OVER (PARTITION BY hadm_id ORDER BY transfertime) AS rn
  FROM `physionet-data.mimiciv_3_1_hosp.services`
)
SELECT
  c.stay_id, c.cancer,
  a.admission_type,
  i.first_careunit,
  p.anchor_year_group,
  s.curr_service AS first_service,
  cs.code_status_raw
FROM coh c
JOIN `physionet-data.mimiciv_3_1_hosp.admissions` a ON a.hadm_id = c.hadm_id
JOIN `physionet-data.mimiciv_3_1_icu.icustays`    i ON i.stay_id = c.stay_id
JOIN `physionet-data.mimiciv_3_1_hosp.patients`   p ON p.subject_id = c.subject_id
LEFT JOIN svc s  ON s.hadm_id  = c.hadm_id  AND s.rn  = 1
LEFT JOIN cs     ON cs.stay_id = c.stay_id AND cs.rn = 1
ORDER BY c.stay_id;

-- Validation: the row count must equal 4,213 + the control cohort size, with no duplicated stay_id
-- SELECT COUNT(*) AS n, COUNT(DISTINCT stay_id) AS distinct_stays FROM (the query above);
