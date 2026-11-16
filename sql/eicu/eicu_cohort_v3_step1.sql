-- eICU-CRD cohort construction, step 1 (Google BigQuery)
-- Replace YOUR_PROJECT with your own BigQuery project id and create the dataset first:
--   CREATE SCHEMA IF NOT EXISTS `YOUR_PROJECT.aki_eicu`;
-- Requires credentialed access to eICU-CRD v2.0 (physionet-data.eicu_crd) and a signed data use agreement.
--
-- This step creates: cancer tiers (A/B/C), substitute covariates and cohort counts.
-- The table contains no outcomes (death, RRT, MAKE), so that the cancer definition can be frozen before outcomes are examined.
--
-- Cancer tiers:
--   A         = ICD-coded cancer (same definition as in the MIMIC-IV analysis)
--   B         = not A, with a neoplasm diagnosis string during the admission, or cancer therapy
--               (chemotherapy or radiotherapy) recorded in the past history
--   C         = not A/B, with only a cancer entry in the past history or an APACHE malignancy flag
--               (may include patients with a history of cancer)
--   ambiguous = only codes for skin/benign neoplasms or a personal history of cancer (excluded)
--   control   = no neoplasm-related signal
--
-- Column names assumed (verify against your copy of eICU-CRD v2.0):
--   patient: unittype, unitadmitsource, hospitaladmitsource
--   apachepatientresult: acutephysiologyscore, apacheversion
--   apachepredvar: metastaticcancer, leukemia, lymphoma, diabetes, cirrhosis, hepaticfailure


CREATE SCHEMA IF NOT EXISTS `YOUR_PROJECT.aki_eicu`;

CREATE OR REPLACE TABLE `YOUR_PROJECT.aki_eicu.cohort_v3` AS
WITH pat AS (
  SELECT
    patientunitstayid, patienthealthsystemstayid, uniquepid, hospitalid, gender,
    unittype, unitadmitsource, hospitaladmitsource, hospitaldischargeyear,
    unitdischargeoffset,
    SAFE_CAST(REPLACE(age, '> ', '') AS INT64) AS age_n,
    ROW_NUMBER() OVER (
      PARTITION BY uniquepid
      ORDER BY hospitaldischargeyear, patienthealthsystemstayid, unitvisitnumber
    ) AS rn
  FROM `physionet-data.eicu_crd.patient`
),
first_icu AS (SELECT * FROM pat WHERE rn = 1 AND age_n >= 18),

-- Creatinine and t0 (creatinine only, KDIGO stage >= 2)
cr AS (
  SELECT patientunitstayid, labresultoffset AS t, labresult AS v
  FROM `physionet-data.eicu_crd.lab`
  WHERE LOWER(labname) = 'creatinine' AND labresult IS NOT NULL
),
cr_base AS (
  SELECT c.patientunitstayid, c.t, c.v, MIN(b.v) AS base
  FROM cr c
  JOIN cr b
    ON b.patientunitstayid = c.patientunitstayid
   AND b.t BETWEEN c.t - 7*24*60 AND c.t
  WHERE c.t >= 0
  GROUP BY c.patientunitstayid, c.t, c.v
),
aki AS (
  SELECT patientunitstayid, MIN(t) AS t0
  FROM cr_base
  WHERE v >= 2 * base
  GROUP BY patientunitstayid
),

-- Diagnosis table: cancer signals, ESRD, comorbidities, sepsis
dx AS (
  SELECT patientunitstayid,
    LOGICAL_OR(
      REGEXP_CONTAINS(IFNULL(icd9code, ''),
        r'(^|,\s*)(1(4\d|5\d|6\d|7[0-2]|7[4-9]|8\d|9\d)|20[0-8])(\.\d+)?\s*(,|$)')
      OR REGEXP_CONTAINS(IFNULL(icd9code, ''),
        r'(^|,\s*)C(0\d|[1-3]\d|4[0-35-9]|[5-7]\d|80|8[1-9]|9[0-6])(\.\w+)?\s*(,|$)')
    ) AS cancer_code,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')),
      r'oncology|cancer|carcinoma|neoplasm|malignan|leukemia|lymphoma|myeloma|sarcoma|metasta')
    ) AS onc_string,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')), r'leukemia|lymphoma|myeloma')
      OR REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)20[0-8](\.\d+)?\s*(,|$)')
      OR REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)C(8[1-9]|9[0-6])(\.\w+)?\s*(,|$)')
    ) AS heme_dx,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')), r'metasta')
      OR REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(19[6-8]|C7[7-9]|C80)')
    ) AS meta_dx,
    LOGICAL_OR(
      REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(1[4-9]\d|20\d|V10)(\.\w+)?\s*(,|$)')
      OR REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(C\d\d|D0\d|D3\d|D4\d|Z85)(\.\w+)?\s*(,|$)')
    ) AS any_neoplasm_code,
    LOGICAL_OR(
      REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(585\.6|V45\.1|N18\.6|Z99\.2)')
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%end stage renal%'
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%chronic dialysis%'
    ) AS esrd_dx,
    -- Comorbidities
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(428|I50)')
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%heart failure%') AS chf_dx,
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(585|N18)')
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%chronic renal%'
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%chronic kidney%') AS ckd_dx,
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(49[0-6]|J4[0-7])')
      OR REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')),
         r'copd|emphysema|chronic obstructive')) AS copd_dx,
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(250|E1[0-3])')
      OR LOWER(IFNULL(diagnosisstring, '')) LIKE '%diabetes%') AS dm_dx,
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(571|572|K7[0-7])')
      OR REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')),
         r'cirrhosis|hepatic failure|hepatic encephalopathy|liver failure')) AS liver_dx,
    LOGICAL_OR(REGEXP_CONTAINS(IFNULL(icd9code, ''), r'(^|,\s*)(038|A4[01]|995\.9[12]|R65\.2)')
      OR REGEXP_CONTAINS(LOWER(IFNULL(diagnosisstring, '')), r'sepsis|septic')) AS sepsis_dx
  FROM `physionet-data.eicu_crd.diagnosis`
  GROUP BY patientunitstayid
),
-- Past history table
ph AS (
  SELECT patientunitstayid,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')),
      r'oncology \(r\)/cancer/')) AS ph_cancer,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')),
      r'oncology \(r\)/cancer therapy/')) AS ph_therapy,
    LOGICAL_OR(LOWER(IFNULL(pasthistorypath, '')) LIKE '%/cancer/metastases/%') AS meta_ph,
    LOGICAL_OR(LOWER(IFNULL(pasthistorypath, '')) LIKE '%/cancer/hematologic malignancy/%') AS heme_ph,
    LOGICAL_OR(LOWER(IFNULL(pasthistorypath, '')) LIKE '%dialysis%') AS esrd_ph,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')), r'heart failure|chf')) AS chf_ph,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')),
      r'renal insufficiency|chronic kidney|chronic renal')) AS ckd_ph,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')),
      r'copd|emphysema|chronic bronchitis')) AS copd_ph,
    LOGICAL_OR(LOWER(IFNULL(pasthistorypath, '')) LIKE '%diabetes%') AS dm_ph,
    LOGICAL_OR(REGEXP_CONTAINS(LOWER(IFNULL(pasthistorypath, '')), r'cirrhosis|hepatic')) AS liver_ph
  FROM `physionet-data.eicu_crd.pasthistory`
  GROUP BY patientunitstayid
),
-- APACHE variables
apv AS (
  SELECT patientunitstayid,
    MAX(IFNULL(metastaticcancer, 0)) AS apv_meta,
    MAX(IFNULL(leukemia, 0))         AS apv_leuk,
    MAX(IFNULL(lymphoma, 0))         AS apv_lymph,
    MAX(IFNULL(diabetes, 0))         AS apv_dm,
    MAX(IFNULL(cirrhosis, 0))        AS apv_cirr,
    MAX(IFNULL(hepaticfailure, 0))   AS apv_hepfail
  FROM `physionet-data.eicu_crd.apachepredvar`
  GROUP BY patientunitstayid
),
aps AS (
  SELECT patientunitstayid,
    MAX(IF(acutephysiologyscore >= 0, acutephysiologyscore, NULL)) AS aps
  FROM `physionet-data.eicu_crd.apachepatientresult`
  WHERE apacheversion = 'IVa'
  GROUP BY patientunitstayid
),
rrt AS (
  SELECT patientunitstayid, MIN(treatmentoffset) AS rrt_t
  FROM `physionet-data.eicu_crd.treatment`
  WHERE LOWER(treatmentstring) LIKE 'renal|dialysis%'
  GROUP BY patientunitstayid
),

base AS (
  SELECT
    f.patientunitstayid, f.hospitalid, f.gender, f.age_n, f.unittype,
    f.unitadmitsource, f.hospitaladmitsource, f.hospitaldischargeyear,
    a.t0,
    IFNULL(d.cancer_code, FALSE) AS cancer_strict,
    IFNULL(d.onc_string, FALSE)  AS onc_string,
    IFNULL(p.ph_therapy, FALSE)  AS ph_therapy,
    (IFNULL(p.ph_cancer, FALSE) OR IFNULL(v.apv_meta, 0) = 1
       OR IFNULL(v.apv_leuk, 0) = 1 OR IFNULL(v.apv_lymph, 0) = 1) AS c_only_signal,
    IFNULL(d.any_neoplasm_code, FALSE) AS any_neoplasm_code,
    (IFNULL(d.heme_dx, FALSE) OR IFNULL(p.heme_ph, FALSE)
       OR IFNULL(v.apv_leuk, 0) = 1 OR IFNULL(v.apv_lymph, 0) = 1) AS heme,
    (IFNULL(d.meta_dx, FALSE) OR IFNULL(p.meta_ph, FALSE) OR IFNULL(v.apv_meta, 0) = 1) AS metastatic,
    (IFNULL(d.esrd_dx, FALSE) OR IFNULL(p.esrd_ph, FALSE)) AS esrd,
    (r.rrt_t IS NOT NULL AND a.t0 IS NOT NULL AND r.rrt_t <= a.t0) AS rrt_before_t0,
    s.aps,
    (IFNULL(d.sepsis_dx, FALSE)) AS sepsis,
    (IFNULL(d.chf_dx, FALSE) OR IFNULL(p.chf_ph, FALSE)) AS chf,
    (IFNULL(d.ckd_dx, FALSE) OR IFNULL(p.ckd_ph, FALSE)) AS ckd,
    (IFNULL(d.copd_dx, FALSE) OR IFNULL(p.copd_ph, FALSE)) AS copd,
    (IFNULL(d.dm_dx, FALSE) OR IFNULL(p.dm_ph, FALSE) OR IFNULL(v.apv_dm, 0) = 1) AS diabetes,
    (IFNULL(d.liver_dx, FALSE) OR IFNULL(p.liver_ph, FALSE)
       OR IFNULL(v.apv_cirr, 0) = 1 OR IFNULL(v.apv_hepfail, 0) = 1) AS liver
  FROM first_icu f
  LEFT JOIN aki a USING (patientunitstayid)
  LEFT JOIN dx d USING (patientunitstayid)
  LEFT JOIN ph p USING (patientunitstayid)
  LEFT JOIN apv v USING (patientunitstayid)
  LEFT JOIN aps s USING (patientunitstayid)
  LEFT JOIN rrt r USING (patientunitstayid)
)
SELECT *,
  CASE
    WHEN cancer_strict THEN 'A'
    WHEN onc_string OR ph_therapy THEN 'B'
    WHEN c_only_signal THEN 'C'
    WHEN any_neoplasm_code THEN 'ambiguous'
    ELSE 'control'
  END AS tier,
  (t0 IS NOT NULL AND NOT esrd AND NOT rrt_before_t0) AS eligible
FROM base;


-- ============================ S1: patients per tier (no outcomes) ============================
SELECT
  tier,
  COUNT(*)                          AS n_eligible,
  COUNTIF(heme)                     AS n_heme,
  COUNTIF(metastatic)               AS n_metastatic,
  COUNT(DISTINCT hospitalid)        AS n_hospitals
FROM `YOUR_PROJECT.aki_eicu.cohort_v3`
WHERE eligible
GROUP BY tier
ORDER BY tier;


-- ============================ S2: covariate availability and prevalence (no outcomes) ============================
SELECT
  tier,
  COUNT(*)                                          AS n_eligible,
  ROUND(AVG(IF(aps IS NOT NULL, 1, 0)), 3)          AS share_with_aps,
  APPROX_QUANTILES(aps, 2)[OFFSET(1)]               AS median_aps,
  APPROX_QUANTILES(age_n, 2)[OFFSET(1)]             AS median_age,
  ROUND(AVG(IF(gender = 'Male', 1, 0)), 3)          AS share_male,
  ROUND(AVG(IF(sepsis, 1, 0)), 3)                   AS share_sepsis,
  ROUND(AVG(IF(chf, 1, 0)), 3)                      AS share_chf,
  ROUND(AVG(IF(ckd, 1, 0)), 3)                      AS share_ckd,
  ROUND(AVG(IF(copd, 1, 0)), 3)                     AS share_copd,
  ROUND(AVG(IF(diabetes, 1, 0)), 3)                 AS share_diabetes,
  ROUND(AVG(IF(liver, 1, 0)), 3)                    AS share_liver,
  ROUND(AVG(IF(unitadmitsource IN ('Operating Room', 'Recovery Room'), 1, 0)), 3) AS share_from_or,
  ROUND(AVG(IF(unitadmitsource = 'Emergency Department', 1, 0)), 3)               AS share_from_ed
FROM `YOUR_PROJECT.aki_eicu.cohort_v3`
WHERE eligible
GROUP BY tier
ORDER BY tier;
