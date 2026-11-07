# Trajectories of acute kidney injury and 90-day mortality in critically ill patients with and without cancer

Code for the manuscript *Trajectories of acute kidney injury and 90-day mortality in critically ill patients with and without cancer: a propensity-matched cohort study* (MIMIC-IV v3.1).

This repository contains code only. **It contains no patient-level data.** MIMIC-IV is a credentialed-access database; its data use agreement does not allow redistribution of the data or of files derived from it (including exports with `stay_id`, `hadm_id` or `subject_id`). Please do not add such files.

## Requirements

- Credentialed access to MIMIC-IV v3.1 on PhysioNet (training completed, data use agreement signed).
- Google BigQuery access to `physionet-data` with the `mimiciv_3_1_hosp`, `mimiciv_3_1_icu` and `mimiciv_3_1_derived` datasets (derived tables used: `kdigo_stages`, `rrt`, `icustay_detail`, `first_day_sofa`, `sepsis3`, `charlson`).
- R 4.6.1 (R 3.6.0 or later is needed to reproduce the fixed-seed results) with the packages `dplyr`, `readr`, `readxl`, `survival`, `MatchIt`, `sandwich`, `lmtest`, `car`, `pROC` and `ggplot2`.
- A BigQuery dataset of your own for the cohort tables. In the SQL files replace `YOUR_PROJECT` by your BigQuery project id; the scripts create and read `YOUR_PROJECT.aki_cancer.cohort` (cancer cohort) and `YOUR_PROJECT.aki_cancer.cohort_ctrl` (non-cancer cohort).

## Folder layout

- `sql/`: S1 to S14, cohort definitions and exports (run in BigQuery)
- `R/`: scripts 05 to 10 and `traj6.R`, the original analysis scripts
- `reproduce/`: `00_export_inputs.R` and `reproduce_all.R`, which regenerate and check every number reported in the manuscript (893 checks)

## Effect measures

Matching uses nine covariates (age, sex, first-day SOFA, sepsis, congestive heart failure, chronic kidney disease, chronic pulmonary disease, diabetes, liver disease). Because admission type, first ICU unit, admitting service and calendar period remained imbalanced after matching, **all outcome models in the manuscript are adjusted for these four variables in addition to the nine matching covariates** (the primary adjustment set, called set B in the code). Sensitivity analyses use the nine matching covariates only (set A) and set B plus the last code status charted at or before t0 (set C); see Section 15 of `reproduce/reproduce_all.R`.

The manuscript reports adjusted **risk ratios (RR)**, estimated by modified Poisson regression with robust variance, and adjusted **risk differences (RD)**, estimated by linear probability models with robust variance (RD in percentage points). In matched analyses standard errors are clustered on matched set; analyses restricted to cancer patients and the unmatched sensitivity analysis use heteroscedasticity-robust standard errors. All of these are computed in `reproduce/reproduce_all.R`. The scripts in `R/` are the earlier analysis scripts. They perform the matching and report odds ratios adjusted for the nine matching covariates only (`07`, `09`, `10` and `traj6.R`); the risk ratios and risk differences in the manuscript, which use the primary adjustment set, are those recomputed by `reproduce/reproduce_all.R`.

## Step 1. SQL (BigQuery)

Run the files in numerical order. Save each exported result as CSV under the name given here, and do not re-save the files in a spreadsheet program.

| SQL | Content | Exported file |
|---|---|---|
| S1 | feasibility counts | (counts only) |
| S2 | cancer cohort (creates `aki_cancer.cohort`) | |
| S3 | non-cancer cohort (creates `aki_cancer.cohort_ctrl`) | |
| S4 | baseline variables, cancer cohort | `baseline.csv` |
| S5 | baseline variables, non-cancer cohort | `baseline_ctrl.csv` |
| S6 | raw creatinine | `creatinine_raw.csv` |
| S7 | KDIGO series (hourly) | `kdigo_series.csv` |
| S8 | daily stage, cancer cohort (used to check the stage derived in R) | `daily_stage_cancer_check.csv` |
| S9 | daily stage, non-cancer cohort | `daily_stage_ctrl.csv` |
| S10 | outcomes, cancer cohort | `outcomes.csv` |
| S11 | outcomes, non-cancer cohort | `outcomes_ctrl.csv` |
| S12 | discharge location (hospice sensitivity analysis) | `adm_hospice.csv` |
| S13 | first-day SOFA and renal SOFA (sensitivity analysis) | `sofa_renal.csv` |
| S14 | extra baseline covariates for both cohorts: admission type, first ICU unit, admitting service, calendar period (anchor-year group) and code status at or before t0 | `extra_confounders.csv` |

Put all exported files in one folder and set `data_dir` to that folder before running the R scripts. `extra_confounders.csv` (required), `adm_hospice.csv` and `sofa_renal.csv` are only read by `reproduce_all.R` (see Step 3). These files contain `stay_id` and must not be committed.

## Step 2. R analysis scripts

Run in one R session, in this order (each script uses objects left by the previous ones):

1. `05_check_data.R` reads the exported files and checks them.
2. `06_rule_based_trajectories.R` derives the daily stage and the rule-based trajectory classes and writes `daily_stage.rds` and `traj_rule_based.rds`.
3. `08_table1_and_figures.R` Table 1 by trajectory, Figure 1 (daily stage) and Figure 2 (Kaplan-Meier from the 72-hour landmark).
4. `09_control_group_comparison.R` adds the non-cancer cohort, performs propensity matching and compares trajectories between cancer and non-cancer patients.
5. `07_followup_checks.R` trajectory analysis in the matched cohort (72-hour landmark, early deaths excluded).
6. `10_cancer_type_subgroups.R` subgroup analysis by cancer type.

`traj6.R` is a stand-alone copy of the six-class trajectory rule (`make_traj6`), which `reproduce_all.R` also defines.

## Step 3. Reproduce and check the reported results

Run in R:

    REPO      <- 'path/to/this/repository'
    PROJ_DIR  <- 'path/to/folder/for/inputs_final.rds'   # also receives reproduce_output/
    ADM_PATH  <- 'path/to/adm_hospice.csv'                # from S12
    SOFA_PATH <- 'path/to/sofa_renal.csv'                 # from S13 (optional)
    EXTRA_PATH <- 'path/to/extra_confounders.csv'         # from S14 (required)
    source(file.path(REPO, 'reproduce', '00_export_inputs.R'))   # run in the session of Step 2: saves inputs_final.rds
    source(file.path(REPO, 'reproduce', 'reproduce_all.R'), echo = FALSE)

`reproduce_all.R` recomputes the tables and compares each number with the value reported in the manuscript; `checks.csv` lists every comparison as OK or DIFF. It covers Tables 1 to 7, Supplementary Figure S1, Supplementary Table S1, the analysis within cancer patients (peak stage and trajectory, likelihood-ratio test, cross-validated AUC), the unmatched analysis, the hospice sensitivity analysis and the sensitivity analysis adjusting for SOFA without its renal component. It also covers the early ICU exit split, the balance of the additional covariates, the unconditional estimates under the three adjustment sets and their E-values (Section 15). On the final data (R 4.6.1) all 893 checks were OK with no DIFF (with `EXTRA_PATH`, `ADM_PATH` and `SOFA_PATH` set; Section 13 not run). `EXTRA_PATH` is required; if `ADM_PATH` or `SOFA_PATH` is missing, Sections 12 and 12b are skipped.

Optional: setting `PE_PATH` to `procedureevents.csv(.gz)` of MIMIC-IV runs Section 13, which looks for the item ids that reproduce the RRT flag. It is not needed for the reported numbers and its output is informational.

The cross-validated AUC uses five folds assigned at random by matched set in the pooled matched sample (`set.seed(1)`), then restricted to cancer patients; the bootstrap interval of the AUC difference uses `set.seed(2)` and 2,000 resamples. The interval depends on this fold assignment.

## Definitions in the code that matter for reading the results

- Index time t0 is the first time with KDIGO stage 2 or higher by creatinine or urine output. RRT is not used to define t0 or any daily stage.
- RRT is the first charted time with active dialysis in the derived `rrt` table, during the ICU stay. The time to RRT (`rrt_h`) is the difference between clock hours (both times truncated to the hour).
- Date of death comes from the post-discharge death records of MIMIC-IV; a patient without a recorded date of death is counted as alive. 90-day death is a date of death within 90 calendar days of the date of t0.
- Discharge to hospice is a `discharge_location` containing "HOSPICE".
- SOFA is the first-day SOFA from the derived `first_day_sofa` table.
- Admission type is grouped as elective or same-day surgical, emergency or urgent, and observation or other; first ICU unit as medical, medical-surgical, cardiac, neuro, and surgical, trauma or other; admitting service as surgical (first service in SURG, CSURG, TSURG, VSURG, NSURG, ORTHO, TRAUM, PSURG, ENT, GU, GYN) or medical/other; calendar period is the three-year `anchor_year_group`. MIMIC-IV has no direct marker of surgery, so these variables are proxies.
- Code status is the last chart entry (item ids 128 and 223758) at or before t0, grouped as full code, DNR/DNI, comfort measures only or not documented; it is documented for only about a third of patients and is used in a sensitivity analysis only.
- E-values follow VanderWeele and Ding (Ann Intern Med 2017;167:268-274).

## Files that must not be committed

`.rds`, `.csv` and `.RData` files, `reproduce_output/` and command histories. The `.gitignore` of this repository already excludes them.

## License and citation

See the `LICENSE` file. Citation: [add the article reference and the Zenodo DOI of the release once available].
