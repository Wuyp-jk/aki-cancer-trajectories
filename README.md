# Trajectories of acute kidney injury and 90-day mortality in critically ill patients with and without cancer

Code for the manuscript *Trajectories of acute kidney injury and 90-day mortality in critically ill patients with and without cancer: a propensity-matched cohort study* (MIMIC-IV v3.1).

This repository contains code only. **It contains no patient-level data.** MIMIC-IV is a credentialed-access database; its data use agreement does not allow redistribution of the data or of files derived from it (including exports with `stay_id`, `hadm_id` or `subject_id`). Please do not add such files.

## Requirements

- Credentialed access to MIMIC-IV v3.1 on PhysioNet (training completed, data use agreement signed).
- Google BigQuery access to `physionet-data` with the `mimiciv_3_1_hosp`, `mimiciv_3_1_icu` and `mimiciv_3_1_derived` datasets (derived tables used: `kdigo_stages`, `rrt`, `icustay_detail`, `first_day_sofa`, `sepsis3`, `charlson`).
- R 4.6.1 with the packages `dplyr`, `readr`, `readxl`, `survival`, `MatchIt`, `sandwich`, `lmtest`, `car`, `pROC` and `ggplot2`.
- A BigQuery dataset of your own for the cohort tables. In the SQL files replace `YOUR_PROJECT` by your BigQuery project id; the scripts create and read `YOUR_PROJECT.aki_cancer.cohort` (cancer cohort) and `YOUR_PROJECT.aki_cancer.cohort_ctrl` (non-cancer cohort).

## Folder layout

- `sql/`: S1 to S11, cohort definitions and exports (run in BigQuery)
- `R/`: scripts 05 to 10 and `traj6.R`, the analysis
- `reproduce/`: `00_export_inputs.R`, `reproduce_all.R` and `section15.R`, which regenerate and check the reported results

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

Put all exported files in one folder and set `data_dir` to that folder before running the R scripts.

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

    REPO     <- 'path/to/this/repository'
    PROJ_DIR <- 'path/to/folder/for/inputs_final.rds'   # also receives reproduce_output/
    ADM_PATH <- 'path/to/file_with_hadm_id_and_discharge_location.csv'
    source(file.path(REPO, 'reproduce', '00_export_inputs.R'))   # run in the session of Step 2: saves inputs_final.rds
    source(file.path(REPO, 'reproduce', 'reproduce_all.R'), echo = FALSE)
    source(file.path(REPO, 'reproduce', 'section15.R'), echo = FALSE)   # analysis within cancer patients

`reproduce_all.R` recomputes the tables and compares each number with the value reported in the manuscript; `checks.csv` lists every comparison as OK or DIFF. The DIFF rows marked as known small differences are kidney outcomes whose definition was rewritten to follow the wording of the Methods; the manuscript reports the recomputed values. `section15.R` should end with 108 OK checks and no DIFF.

## Definitions in the code that matter for reading the results

- Index time t0 is the first time with KDIGO stage 2 or higher by creatinine or urine output. RRT is not used to define t0 or any daily stage.
- RRT is the first charted time with active dialysis in the derived `rrt` table, during the ICU stay.
- Date of death comes from the post-discharge death records of MIMIC-IV; a patient without a recorded date of death is counted as alive. 90-day death is a date of death within 90 calendar days of the date of t0.

## Files that must not be committed

`.rds`, `.csv` and `.RData` files, `reproduce_output/` and command histories. The `.gitignore` of this repository already excludes them.

## License and citation

Released under the MIT License (see LICENSE)
