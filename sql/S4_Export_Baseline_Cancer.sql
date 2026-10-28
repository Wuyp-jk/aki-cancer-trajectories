-- ============================================================
-- 04b  Baseline Table (Cancer Cohort, One row per patient)
-- Prerequisite: 04a has been executed
-- Export: Click "Save Results" -> CSV, name the file baseline.csv
-- ============================================================

WITH ncr AS (
  -- Number of creatinine measurements per patient within 7 days after t0 (during ICU stay), same criteria as before
  SELECT c.stay_id, COUNT(*) AS n_cr
  FROM `YOUR_PROJECT.aki_cancer.cohort` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON l.hadm_id = c.hadm_id
   AND l.itemid = 50912
   AND l.valuenum IS NOT NULL
   AND l.charttime BETWEEN c.t0
       AND LEAST(DATETIME_ADD(c.t0, INTERVAL 7 DAY), c.icu_outtime)
  GROUP BY c.stay_id
),
-- Baseline creatinine: Minimum creatinine within 7 days before t0 (including t0); null if not measured during this period
bl AS (
  SELECT c.stay_id, MIN(l.valuenum) AS creat_baseline_7d
  FROM `YOUR_PROJECT.aki_cancer.cohort` c
  JOIN `physionet-data.mimiciv_3_1_hosp.labevents` l
    ON l.hadm_id = c.hadm_id
   AND l.itemid = 50912
   AND l.valuenum IS NOT NULL
   AND l.charttime BETWEEN DATETIME_SUB(c.t0, INTERVAL 7 DAY) AND c.t0
  GROUP BY c.stay_id
)
SELECT
  c.stay_id, c.subject_id, c.hadm_id,
  c.gender, c.race, c.admission_age,
  c.solid, c.heme,
  COALESCE(n.n_cr, 0)            AS n_cr_7d,
  (COALESCE(n.n_cr, 0) >= 3)     AS ge3_cr,        
  b.creat_baseline_7d,
  sf.sofa                        AS sofa_day1,
  COALESCE(s3.sepsis3, FALSE)    AS sepsis3,
  -- Charlson comorbidities (including metastatic_solid_tumor, renal_disease, etc.)
  ch.* EXCEPT (subject_id, hadm_id)
FROM `YOUR_PROJECT.aki_cancer.cohort` c
LEFT JOIN ncr n  ON n.stay_id = c.stay_id
LEFT JOIN bl  b  ON b.stay_id = c.stay_id
LEFT JOIN `physionet-data.mimiciv_3_1_derived.first_day_sofa` sf ON sf.stay_id = c.stay_id
LEFT JOIN `physionet-data.mimiciv_3_1_derived.sepsis3`        s3 ON s3.stay_id = c.stay_id
LEFT JOIN `physionet-data.mimiciv_3_1_derived.charlson`       ch ON ch.hadm_id = c.hadm_id
ORDER BY c.stay_id;
