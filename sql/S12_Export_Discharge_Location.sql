-- S12: discharge location for every admission (used for the hospice sensitivity analysis).
-- Reads the public MIMIC-IV table only; no cohort table is needed.
-- Save the result as adm_hospice.csv (columns: hadm_id, discharge_location).
SELECT hadm_id, discharge_location
FROM `physionet-data.mimiciv_3_1_hosp.admissions`
