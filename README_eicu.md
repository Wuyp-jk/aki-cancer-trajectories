## External replication in eICU-CRD

This part of the repository reproduces the external replication reported in the manuscript
(supplementary Methods, Table 9 and Tables S6-S10). The replication uses the eICU Collaborative
Research Database (eICU-CRD v2.0). Patient-level data cannot be shared; the scripts below rebuild the
analysis files for researchers with credentialed access, and `results/eicu/` holds the aggregate results
reported in the paper.

### Requirements
- Credentialed access to eICU-CRD v2.0 on PhysioNet, with the data use agreement signed, and access to
  `physionet-data.eicu_crd` on Google BigQuery.
- A Google BigQuery project of your own (replace `YOUR_PROJECT` in the SQL files) and a dataset named `aki_eicu`.
- R >= 4.0 with the packages `MatchIt`, `sandwich`, `lmtest` and `dplyr`.

### Run order
1. `sql/eicu/eicu_cohort_v3_step1.sql` creates `YOUR_PROJECT.aki_eicu.cohort_v3` (cancer tiers A/B/C, substitute
   covariates; no outcomes) and prints the number of patients per tier.
2. `sql/eicu/eicu_cohort_v3_step2_fixed.sql` creates `YOUR_PROJECT.aki_eicu.analysis_v3` (daily KDIGO stage,
   trajectory class and outcomes).
3. Export `analysis_v3` (`SELECT * FROM YOUR_PROJECT.aki_eicu.analysis_v3`, save the result as CSV). Keep the file
   on your own computer; do not commit or share it.
4. In R, set `DATA_PATH` to the exported CSV (and optionally `OUT_DIR`) and run `R/eicu/eicu_analysis.R`. The script
   writes seven aggregate CSV files (counts below 10 are shown as "<10").

### Cohort definition and differences from the MIMIC-IV analysis
- Adults in their first ICU stay of the first hospitalization with creatinine-based KDIGO stage >= 2 AKI; end-stage
  renal disease, chronic dialysis and renal replacement therapy before t0 excluded. The exclusion for an urgent
  dialysis indication before t0 is not applied.
- AKI is staged by serum creatinine only (urine output is charted too sparsely).
- Only in-hospital outcomes are available (no 90-day mortality). The class without creatinine on days 1-2 is named
  `no_data` and is not equivalent to the "early ICU exit" class of the MIMIC-IV analysis.
- Cancer is identified by a tiered algorithm: tier A (ICD-coded cancer), tier B (neoplasm diagnosis during the
  admission or cancer therapy in the past history) and tier C (past-history cancer entry or APACHE malignancy flag only).
  The primary analysis uses tiers A to C; tiers A+B and tier A alone are sensitivity analyses.
- The APACHE IV acute physiology score replaces SOFA; sepsis and comorbidities are derived from diagnosis, past
  history and APACHE variables.
- `make_hosp`, `make_hosp_strict` and `nonrecovery_raw` in `analysis_v3` were not analyzed (the control-group
  frequency of MAKE-hosp, about 54%, was implausible).

### Results files (`results/eicu/`)
| File | Content |
|---|---|
| `01_estimates_all.csv` | RR and RD of cancer vs control for all analyses and cancer definitions |
| `02_comparison_with_mimic.csv` | eICU-CRD estimates next to the MIMIC-IV estimates (Table 9) |
| `03_interactions.csv` | Interaction tests (cancer x trajectory, cancer x cancer type) |
| `04_gradient_within_cancer.csv` | In-hospital death by trajectory within cancer patients (Table S9) |
| `05_trajectory_distribution.csv`, `05b_...` | Trajectory distribution and chi-square test (Table S8) |
| `06_balance.csv` | Standardized mean differences before and after matching (Table S10) |

Note: in `02_comparison_with_mimic.csv` the row `all_matched_rrt_any` is not a like-for-like comparison, because the MIMIC-IV RRT estimate comes from the trajectory sample; Table 9 uses the row `traj_sample_rrt_any` of `01_estimates_all.csv` instead. The comparison with MIMIC-IV uses the reference values transcribed from
Tables 5, 6 and 8 of the manuscript, and the z-test P values reported in Table 9 are computed from the confidence
intervals (post hoc, not adjusted for multiple testing).
