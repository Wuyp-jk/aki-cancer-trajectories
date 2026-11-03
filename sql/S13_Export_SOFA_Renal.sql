-- S13: first-day SOFA and its renal component (used for the sensitivity analysis that
-- adjusts for SOFA without the renal term, worsening vs severe-from-start AKI).
-- Reads the public MIMIC-IV derived table only; no cohort table is needed.
-- Save the result as sofa_renal.csv (columns: stay_id, sofa, renal).
SELECT stay_id, sofa, renal
FROM `physionet-data.mimiciv_3_1_derived.first_day_sofa`
