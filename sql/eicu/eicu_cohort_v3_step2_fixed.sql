-- eICU-CRD cohort construction, step 2 (Google BigQuery)
-- Prerequisite: step 1 has been run and the table YOUR_PROJECT.aki_eicu.cohort_v3 exists.
-- Creates the table analysis_v3: daily KDIGO stage (creatinine only), trajectory class
-- (rules and order of Table 1 of the manuscript) and outcomes.
--
-- Differences from the MIMIC-IV analysis:
--   1. AKI is staged by creatinine only (urine output is charted too sparsely in eICU-CRD).
--   2. eICU-CRD has no post-discharge death dates, so the outcome is in-hospital death; there is no 90-day mortality.
--   3. A patient without a creatinine measurement on days 1-2 may still be in the ICU; this class is named no_data
--      (it is not equivalent to the "early ICU exit" class of the MIMIC-IV analysis).
--   4. Baseline creatinine is the lowest value in the 7 days up to and including t0 (as in the MIMIC-IV analysis);
--      eICU-CRD offsets are in minutes.
--   5. The exclusion for an urgent dialysis indication before t0 is not applied.
--   6. The creatinine-based columns make_hosp, make_hosp_strict and nonrecovery_raw are kept for transparency but were
--      NOT analyzed in the manuscript, because the control-group frequency (about 54% for MAKE-hosp) was implausible.
--
-- Column names assumed (verify against your copy): patient.hospitaldischargelocation, patient.unitdischargestatus


CREATE OR REPLACE TABLE `YOUR_PROJECT.aki_eicu.analysis_v3` AS
WITH c AS (
  SELECT * FROM `YOUR_PROJECT.aki_eicu.cohort_v3`
  WHERE eligible AND tier IN ('A', 'B', 'C', 'control')
),
pat AS (
  SELECT patientunitstayid, hospitaldischargestatus, hospitaldischargeoffset,
         hospitaldischargelocation, unitdischargeoffset
  FROM `physionet-data.eicu_crd.patient`
),
cr AS (
  SELECT patientunitstayid, labresultoffset AS t, labresult AS v
  FROM `physionet-data.eicu_crd.lab`
  WHERE LOWER(labname) = 'creatinine' AND labresult IS NOT NULL
),
base AS (
  SELECT c.patientunitstayid, MIN(r.v) AS base_cr
  FROM c
  JOIN cr r ON r.patientunitstayid = c.patientunitstayid
           AND r.t BETWEEN c.t0 - 7*24*60 AND c.t0
  GROUP BY c.patientunitstayid
),
meas AS (
  SELECT c.patientunitstayid,
         CAST(FLOOR((r.t - c.t0) / 1440) AS INT64) AS day,
         IF(r.v >= 3 * b.base_cr OR r.v >= 4.0, 3, IF(r.v >= 2 * b.base_cr, 2, 1)) AS stg
  FROM c
  JOIN base b USING (patientunitstayid)
  JOIN cr r ON r.patientunitstayid = c.patientunitstayid
  WHERE r.t >= c.t0 AND r.t < c.t0 + 3*1440
),
daily AS (
  SELECT patientunitstayid,
         MAX(IF(day = 0, stg, NULL)) AS d0,
         MAX(IF(day = 1, stg, NULL)) AS d1,
         MAX(IF(day = 2, stg, NULL)) AS d2
  FROM meas
  GROUP BY patientunitstayid
),
last_cr AS (
  SELECT patientunitstayid,
         ARRAY_AGG(STRUCT(t, v) ORDER BY t DESC LIMIT 1)[OFFSET(0)] AS l
  FROM (
    SELECT r.patientunitstayid, r.t, r.v
    FROM cr r
    JOIN pat p USING (patientunitstayid)
    WHERE r.t <= p.hospitaldischargeoffset
  )
  GROUP BY patientunitstayid
),
rrt AS (
  SELECT patientunitstayid, MIN(treatmentoffset) AS rrt_t
  FROM `physionet-data.eicu_crd.treatment`
  WHERE LOWER(treatmentstring) LIKE 'renal|dialysis%'
  GROUP BY patientunitstayid
),
full_tbl AS (
  SELECT
    c.patientunitstayid, c.hospitalid, c.tier, (c.tier IN ('A', 'B', 'C')) AS cancer,
    c.heme, c.metastatic, c.age_n, c.gender, c.aps, c.sepsis, c.chf, c.ckd, c.copd,
    c.diabetes, c.liver, c.unittype, c.unitadmitsource, c.t0,
    b.base_cr, d.d0, d.d1, d.d2,
    p.hospitaldischargeoffset, p.unitdischargeoffset, p.hospitaldischargelocation,
    (p.hospitaldischargestatus = 'Expired') AS hosp_death,
    (p.hospitaldischargestatus = 'Expired'
       AND p.hospitaldischargeoffset <= c.t0 + 72*60) AS early_death,
    r.rrt_t,
    (r.rrt_t IS NOT NULL) AS rrt_any,
    (r.rrt_t IS NOT NULL AND r.rrt_t <= c.t0 + 72*60) AS rrt_first72h,
    (r.rrt_t IS NOT NULL AND r.rrt_t > c.t0 + 72*60) AS rrt_after72h,
    lc.l.v AS last_cr, lc.l.t AS last_cr_t
  FROM c
  JOIN base b USING (patientunitstayid)
  JOIN pat p USING (patientunitstayid)
  LEFT JOIN daily d USING (patientunitstayid)
  LEFT JOIN rrt r USING (patientunitstayid)
  LEFT JOIN last_cr lc USING (patientunitstayid)
)
SELECT *,
  -- Trajectory class: rules of Table 1 applied in order; the first rule met decides the class
  CASE
    WHEN early_death THEN '1_early_death'
    WHEN d1 IS NULL AND d2 IS NULL THEN '2_no_data'
    WHEN IFNULL(d1, 0) < 2 AND IFNULL(d2, 0) < 2 THEN '3_transient'
    WHEN d0 = 3 THEN '6_severe_from_start'
    WHEN IFNULL(d1, 0) = 3 OR IFNULL(d2, 0) = 3 THEN '5_worsening'
    ELSE '4_persistent_stage2'
  END AS traj,
  -- Kidney outcomes (in hospital)
  -- The discharge creatinine must be measured within 72 hours before discharge; otherwise the creatinine at t0 would be mistaken for it
  (hosp_death OR rrt_any OR (last_cr IS NOT NULL
     AND last_cr_t >= hospitaldischargeoffset - 72*60
     AND last_cr >= 2 * base_cr)) AS make_hosp,
  -- Strict version: the creatinine must also be measured after t0
  (hosp_death OR rrt_any OR (last_cr IS NOT NULL AND last_cr_t > t0
     AND last_cr_t >= hospitaldischargeoffset - 72*60
     AND last_cr >= 2 * base_cr)) AS make_hosp_strict,
  (NOT hosp_death AND base_cr >= 0.2 AND last_cr IS NOT NULL
     AND last_cr_t >= hospitaldischargeoffset - 72*60) AS nr_evaluable,
  (rrt_any OR (last_cr IS NOT NULL AND last_cr > 1.5 * base_cr)) AS nonrecovery_raw,
  (NOT hosp_death AND base_cr >= 0.2 AND last_cr IS NOT NULL AND last_cr_t > t0
     AND last_cr_t >= hospitaldischargeoffset - 72*60) AS nr_evaluable_strict
FROM full_tbl;


-- ============================ S3: trajectory distribution and crude in-hospital outcomes, cancer vs control ============================
SELECT
  cancer, traj,
  COUNT(*)                                    AS n,
  ROUND(AVG(IF(hosp_death, 1, 0)), 3)         AS hosp_death_rate,
  ROUND(AVG(IF(make_hosp, 1, 0)), 3)          AS make_hosp_rate,
  ROUND(AVG(IF(make_hosp_strict, 1, 0)), 3)   AS make_hosp_strict_rate,
  ROUND(AVG(IF(rrt_any, 1, 0)), 3)            AS rrt_rate
FROM `YOUR_PROJECT.aki_eicu.analysis_v3`
GROUP BY cancer, traj
ORDER BY cancer DESC, traj;


-- ============================ S4: crude outcomes overall and by cancer tier ============================
SELECT
  tier,
  COUNT(*)                                    AS n,
  ROUND(AVG(IF(hosp_death, 1, 0)), 3)         AS hosp_death_rate,
  ROUND(AVG(IF(early_death, 1, 0)), 3)        AS early_death_rate,
  ROUND(AVG(IF(make_hosp, 1, 0)), 3)          AS make_hosp_rate,
  ROUND(AVG(IF(make_hosp_strict, 1, 0)), 3)   AS make_hosp_strict_rate,
  ROUND(AVG(IF(rrt_any, 1, 0)), 3)            AS rrt_rate,
  COUNTIF(nr_evaluable)                       AS n_nr_evaluable,
  COUNTIF(nr_evaluable_strict)                AS n_nr_evaluable_strict
FROM `YOUR_PROJECT.aki_eicu.analysis_v3`
GROUP BY tier
ORDER BY tier;
