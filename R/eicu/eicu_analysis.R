# =====================================================================
# External replication in eICU-CRD (run locally)
#
# Input : the table analysis_v3 exported from BigQuery (patient-level data; keep it on your own
#         computer and do not share it: the eICU-CRD data use agreement does not allow sharing
#         patient-level data).
# Output: aggregate CSV files in OUT_DIR (counts, RR, RD, confidence intervals, P values;
#         counts below 10 are shown as "<10").
#
# Primary analysis    : cancer = tiers A+B+C vs controls; 1:2 nearest-neighbor matching on the
#                       logit propensity score (caliper 0.2 SD).
# Sensitivity analyses: cancer = tiers A+B; tier A only.
# Outcomes            : in-hospital death (primary) and renal replacement therapy (secondary).
#                       Creatinine-based MAKE-hosp and non-recovery are not analyzed.
#
# Requires R >= 4.0 and the packages MatchIt, sandwich, lmtest and dplyr.
# Usage: set DATA_PATH (and optionally OUT_DIR) in the console, or edit the defaults below, then
#        run the whole script.
# =====================================================================

suppressPackageStartupMessages({
  library(MatchIt); library(sandwich); library(lmtest); library(dplyr)
})

# Values already set in the console (e.g. DATA_PATH <- file.choose()) are not overwritten here
if (!exists("DATA_PATH")) DATA_PATH <- "analysis_v3.csv"
if (!exists("OUT_DIR"))   OUT_DIR   <- "eicu_results"
dir.create(OUT_DIR, showWarnings = FALSE)
set.seed(1)

# ---------------------------- Read and prepare data ----------------------------
tobool <- function(x) if (is.logical(x)) x else tolower(as.character(x)) %in% c("true", "t", "1")
if (!file.exists(DATA_PATH)) {
  stop("File not found: ", DATA_PATH, "\nCurrent working directory: ", getwd(),
       "\nSet DATA_PATH to the full path of the exported CSV, e.g. \"/path/to/analysis_v3.csv\"")
}
d <- read.csv(DATA_PATH, stringsAsFactors = FALSE)
need <- c("tier", "cancer", "heme", "metastatic", "age_n", "gender", "aps", "sepsis", "chf", "ckd",
          "copd", "diabetes", "liver", "unittype", "unitadmitsource", "traj", "hosp_death",
          "early_death", "rrt_any", "hospitalid")
if (length(setdiff(need, names(d))) > 0) {
  stop("The CSV is missing these columns: ", paste(setdiff(need, names(d)), collapse = ", "),
       "\nExport all columns of analysis_v3 (SELECT *).")
}
cat("Read", nrow(d), "rows and", ncol(d), "columns\n")

for (v in c("cancer", "heme", "metastatic", "sepsis", "chf", "ckd", "copd", "diabetes", "liver",
            "hosp_death", "early_death", "rrt_any", "rrt_first72h", "rrt_after72h")) {
  d[[v]] <- tobool(d[[v]])
}
d$aps <- suppressWarnings(as.numeric(d$aps))

d <- d %>% mutate(
  age         = age_n,
  male        = as.integer(gender == "Male"),
  aps_missing = as.integer(is.na(aps)),
  aps_f       = ifelse(is.na(aps), median(aps, na.rm = TRUE), aps),
  sepsis = as.integer(sepsis), chf = as.integer(chf), ckd = as.integer(ckd),
  copd = as.integer(copd), diabetes = as.integer(diabetes), liver = as.integer(liver),
  icu_unit = case_when(
    unittype %in% c("Cardiac ICU", "CCU-CTICU", "CTICU", "CSICU") ~ "cardiac",
    unittype == "Med-Surg ICU" ~ "medsurg",
    unittype == "MICU"         ~ "medical",
    unittype == "Neuro ICU"    ~ "neuro",
    TRUE                       ~ "surgical_other"),
  adm_source = case_when(
    unitadmitsource == "Emergency Department"               ~ "ED",
    unitadmitsource %in% c("Operating Room", "Recovery Room") ~ "OR",
    TRUE                                                    ~ "other"),
  ctype = case_when(
    !cancer    ~ NA_character_,
    heme       ~ "hematologic",
    metastatic ~ "metastatic_solid",
    TRUE       ~ "nonmetastatic_solid")
)

TRAJ_KEEP <- c("3_transient", "4_persistent_stage2", "5_worsening", "6_severe_from_start")

# ---------------------------- Helper functions ----------------------------
fmt_n <- function(n) ifelse(n > 0 & n < 10, "<10", as.character(n))

# Cancer vs control: RR (modified Poisson) and RD (linear probability model, percentage points)
est_effect <- function(dat, outcome, rhs, cluster = NULL) {
  dat$y <- as.integer(dat[[outcome]])
  f <- as.formula(paste("y ~ cancer_i +", rhs))
  vc <- function(m) if (is.null(cluster)) vcovHC(m, type = "HC1") else vcovCL(m, cluster = cluster, type = "HC1")
  out <- tryCatch({
    fp <- glm(f, data = dat, family = poisson(link = "log"))
    cp <- coeftest(fp, vcov = vc(fp))["cancer_i", ]
    fl <- lm(f, data = dat)
    cl <- coeftest(fl, vcov = vc(fl))["cancer_i", ]
    data.frame(
      n_cancer = sum(dat$cancer_i == 1), n_ctrl = sum(dat$cancer_i == 0),
      risk_cancer = round(100 * mean(dat$y[dat$cancer_i == 1]), 1),
      risk_ctrl   = round(100 * mean(dat$y[dat$cancer_i == 0]), 1),
      rr = exp(cp[1]), rr_lo = exp(cp[1] - 1.96 * cp[2]), rr_hi = exp(cp[1] + 1.96 * cp[2]),
      rd = 100 * cl[1], rd_lo = 100 * (cl[1] - 1.96 * cl[2]), rd_hi = 100 * (cl[1] + 1.96 * cl[2]))
  }, error = function(e) {
    data.frame(n_cancer = sum(dat$cancer_i == 1), n_ctrl = sum(dat$cancer_i == 0),
               risk_cancer = NA, risk_ctrl = NA, rr = NA, rr_lo = NA, rr_hi = NA,
               rd = NA, rd_lo = NA, rd_hi = NA)
  })
  out
}

# Interaction tests on the multiplicative and additive scales (cancer x stratifying variable)
inter_test <- function(dat, outcome, strat, rhs) {
  dat$y <- as.integer(dat[[outcome]]); dat$s <- factor(dat[[strat]])
  f0 <- as.formula(paste("y ~ cancer_i + s +", rhs))
  f1 <- as.formula(paste("y ~ cancer_i * s +", rhs))
  vcf <- function(x) vcovCL(x, cluster = dat$ms, type = "HC1")
  tryCatch({
    p0 <- glm(f0, data = dat, family = poisson); p1 <- glm(f1, data = dat, family = poisson)
    l0 <- lm(f0, data = dat);                     l1 <- lm(f1, data = dat)
    c(p_mult = waldtest(p0, p1, vcov = vcf, test = "Chisq")[["Pr(>Chisq)"]][2],
      p_add  = waldtest(l0, l1, vcov = vcf, test = "Chisq")[["Pr(>Chisq)"]][2])
  }, error = function(e) c(p_mult = NA, p_add = NA))
}

# Within cancer patients: RR / RD of each trajectory vs the reference class (HC1 robust standard errors)
est_levels <- function(dat, outcome, var, ref, rhs) {
  dat$y <- as.integer(dat[[outcome]]); dat$g <- relevel(factor(dat[[var]]), ref = ref)
  f <- as.formula(paste("y ~ g +", rhs))
  tryCatch({
    fp <- glm(f, data = dat, family = poisson); cp <- coeftest(fp, vcov = vcovHC(fp, type = "HC1"))
    fl <- lm(f, data = dat);                    cl <- coeftest(fl, vcov = vcovHC(fl, type = "HC1"))
    rows <- grep("^g", rownames(cp))
    data.frame(level = sub("^g", "", rownames(cp)[rows]),
      n = as.vector(table(dat$g)[sub("^g", "", rownames(cp)[rows])]),
      risk = round(100 * tapply(dat$y, dat$g, mean)[sub("^g", "", rownames(cp)[rows])], 1),
      rr = exp(cp[rows, 1]), rr_lo = exp(cp[rows, 1] - 1.96 * cp[rows, 2]),
      rr_hi = exp(cp[rows, 1] + 1.96 * cp[rows, 2]),
      rd = 100 * cl[rows, 1], rd_lo = 100 * (cl[rows, 1] - 1.96 * cl[rows, 2]),
      rd_hi = 100 * (cl[rows, 1] + 1.96 * cl[rows, 2]), row.names = NULL)
  }, error = function(e) NULL)
}

smd <- function(x, g) {
  s <- sqrt((var(x[g == 1]) + var(x[g == 0])) / 2)
  if (is.na(s) || s == 0) return(0)
  (mean(x[g == 1]) - mean(x[g == 0])) / s
}

# ---------------------------- Matching ----------------------------
do_match <- function(dat) {
  mv <- c("age", "male", "aps_f", "sepsis", "chf", "ckd", "copd", "diabetes", "liver")
  if (any(dat$aps_missing == 1)) mv <- c(mv, "aps_missing")
  m <- matchit(as.formula(paste("cancer_i ~", paste(mv, collapse = " + "))),
               data = dat, method = "nearest", distance = "glm", link = "linear.logit",
               ratio = 2, caliper = 0.2, std.caliper = TRUE, replace = FALSE)
  md <- match.data(m)
  md$ms <- as.character(md$subclass)
  ct <- tapply(ifelse(md$cancer_i == 1, md$ctype, NA_character_), md$ms,
               function(x) x[!is.na(x)][1])
  md$ctype_m <- unname(ct[md$ms])
  list(dat = dat, md = md, mv = mv)
}

RHS <- "age + male + aps_f + aps_missing + sepsis + chf + ckd + copd + diabetes + liver + icu_unit + adm_source"

# ---------------------------- Main analysis ----------------------------
analyze <- function(tiers, tag, full = TRUE) {
  dat <- d[d$tier %in% c(tiers, "control"), ]
  dat$cancer_i <- as.integer(dat$tier != "control")
  rhs <- if (any(dat$aps_missing == 1)) RHS else sub("aps_missing \\+ ", "", RHS)
  mm <- do_match(dat); md <- mm$md
  cat(sprintf("[%s] cancer %d, control pool %d; after matching: cancer %d, controls %d\n", tag,
              sum(dat$cancer_i == 1), sum(dat$cancer_i == 0),
              sum(md$cancer_i == 1), sum(md$cancer_i == 0)))
  res <- list()
  mdt <- md[md$traj %in% TRAJ_KEEP, ]

  add <- function(name, df) { df$analysis <- name; df$cohort <- tag; res[[length(res) + 1]] <<- df }

  # 1) All matched patients (including early deaths): main estimates
  add("all_matched_hosp_death", est_effect(md, "hosp_death", rhs, md$ms))
  add("all_matched_rrt_any",    est_effect(md, "rrt_any", rhs, md$ms))
  # 2) Trajectory sample overall (adjusted for trajectory)
  add("traj_sample_hosp_death", est_effect(mdt, "hosp_death", paste(rhs, "+ traj"), mdt$ms))
  add("traj_sample_rrt_any",    est_effect(mdt, "rrt_any",   paste(rhs, "+ traj"), mdt$ms))
  # 3) Sensitivity analysis clustering on hospital (matched sets may not capture hospital-level correlation)
  add("all_matched_hosp_death_hospital_cluster", est_effect(md, "hosp_death", rhs, md$hospitalid))

  out <- list(res = res, md = md, mdt = mdt, mm = mm, rhs = rhs)
  if (!full) return(out)

  # 4) Within each trajectory
  for (tj in TRAJ_KEEP) {
    s <- mdt[mdt$traj == tj, ]
    add(paste0("traj_", tj, "_hosp_death"), est_effect(s, "hosp_death", rhs, s$ms))
  }
  # 5) By cancer type
  for (ty in c("hematologic", "nonmetastatic_solid", "metastatic_solid")) {
    s <- mdt[!is.na(mdt$ctype_m) & mdt$ctype_m == ty, ]
    add(paste0("type_", ty, "_hosp_death"), est_effect(s, "hosp_death", rhs, s$ms))
  }
  list(res = res, md = md, mdt = mdt, mm = mm, rhs = rhs)
}

primary <- analyze(c("A", "B", "C"), "ABC", full = TRUE)
sens_AB <- analyze(c("A", "B"), "AB", full = FALSE)
sens_A  <- analyze("A", "A", full = FALSE)

md <- primary$md; mdt <- primary$mdt; rhs <- primary$rhs

# ---------------------------- Interaction tests ----------------------------
it_traj <- inter_test(mdt, "hosp_death", "traj", rhs)
mdt_type <- mdt[!is.na(mdt$ctype_m), ]
it_type <- inter_test(mdt_type, "hosp_death", "ctype_m", rhs)
interactions <- data.frame(
  test = c("cancer x trajectory", "cancer x cancer type"),
  p_multiplicative = c(it_traj["p_mult"], it_type["p_mult"]),
  p_additive       = c(it_traj["p_add"],  it_type["p_add"]))

# ---------------------------- Gradient across trajectories within cancer patients ----------------------------
canc <- mdt[mdt$cancer_i == 1, ]
gradient <- est_levels(canc, "hosp_death", "traj", "3_transient", rhs)

# ---------------------------- Trajectory distribution (descriptive) ----------------------------
tab <- table(md$traj, ifelse(md$cancer_i == 1, "cancer", "control"))
dist <- data.frame(traj = rownames(tab),
                   cancer_n = fmt_n(tab[, "cancer"]),
                   cancer_pct = round(100 * tab[, "cancer"] / sum(tab[, "cancer"]), 1),
                   control_n = fmt_n(tab[, "control"]),
                   control_pct = round(100 * tab[, "control"] / sum(tab[, "control"]), 1))
chi <- suppressWarnings(chisq.test(tab))
dist_p <- data.frame(chisq = unname(chi$statistic), df = unname(chi$parameter), p = chi$p.value)

# ---------------------------- Covariate balance ----------------------------
bal_vars <- c("age", "male", "aps_f", "sepsis", "chf", "ckd", "copd", "diabetes", "liver")
all_dat <- primary$mm$dat
balance <- data.frame(
  variable = bal_vars,
  smd_before = sapply(bal_vars, function(v) round(smd(all_dat[[v]], all_dat$cancer_i), 3)),
  smd_after  = sapply(bal_vars, function(v) round(smd(md[[v]], md$cancer_i), 3)), row.names = NULL)

# ---------------------------- Collect all estimates ----------------------------
est_all <- do.call(rbind, c(primary$res, sens_AB$res, sens_A$res))
num_cols <- c("rr", "rr_lo", "rr_hi", "rd", "rd_lo", "rd_hi")
est_all[num_cols] <- lapply(est_all[num_cols], function(x) round(x, 2))
est_all$n_cancer <- fmt_n(est_all$n_cancer); est_all$n_ctrl <- fmt_n(est_all$n_ctrl)
est_all <- est_all[, c("cohort", "analysis", "n_cancer", "n_ctrl", "risk_cancer", "risk_ctrl",
                       "rr", "rr_lo", "rr_hi", "rd", "rd_lo", "rd_hi")]

# ---------------------------- Comparison with MIMIC-IV ----------------------------
# MIMIC-IV reference estimates transcribed from Tables 5, 6 and 8 of the manuscript
ref <- data.frame(
  analysis = c("all_matched_hosp_death", "traj_sample_hosp_death",
               "traj_3_transient_hosp_death", "traj_4_persistent_stage2_hosp_death",
               "traj_5_worsening_hosp_death", "traj_6_severe_from_start_hosp_death",
               "type_hematologic_hosp_death", "type_nonmetastatic_solid_hosp_death",
               "type_metastatic_solid_hosp_death", "all_matched_rrt_any"),
  m_rr    = c(1.51, 1.39, 1.74, 1.59, 1.02, 1.22, 1.61, 0.98, 1.74, 0.94),
  m_rr_lo = c(1.40, 1.26, 1.39, 1.32, 0.81, 1.03, 1.26, 0.82, 1.50, 0.78),
  m_rr_hi = c(1.63, 1.54, 2.18, 1.91, 1.30, 1.45, 2.05, 1.18, 2.03, 1.13),
  m_rd    = c(7.6, 5.4, 5.0, 7.2, 1.4, 5.1, 8.6, -0.4, 9.6, -0.7),
  m_rd_lo = c(6.1, 3.6, 2.7, 4.1, -6.2, 0.3, 4.2, -3.1, 6.7, -1.7),
  m_rd_hi = c(9.0, 7.2, 7.3, 10.2, 9.0, 9.8, 13.0, 2.3, 12.4, 0.4))
# Note: the MIMIC-IV RRT estimate comes from the trajectory sample (Table 8 of the manuscript); the row all_matched_rrt_any is shown for direction only.

cmp <- merge(est_all[est_all$cohort == "ABC", ], ref, by = "analysis")
cmp$same_side_of_null_rr <- (cmp$rr > 1) == (cmp$m_rr > 1)
cmp$rr_in_mimic_ci <- cmp$rr >= cmp$m_rr_lo & cmp$rr <= cmp$m_rr_hi
cmp$rd_in_mimic_ci <- cmp$rd >= cmp$m_rd_lo & cmp$rd <= cmp$m_rd_hi
cmp$ci_overlap_rd  <- cmp$rd_lo <= cmp$m_rd_hi & cmp$rd_hi >= cmp$m_rd_lo
cmp <- cmp[, c("analysis", "rr", "rr_lo", "rr_hi", "m_rr", "m_rr_lo", "m_rr_hi",
               "rd", "rd_lo", "rd_hi", "m_rd", "m_rd_lo", "m_rd_hi",
               "same_side_of_null_rr", "rr_in_mimic_ci", "rd_in_mimic_ci", "ci_overlap_rd")]

# ---------------------------- Write output (aggregate results only) ----------------------------
write.csv(est_all,      file.path(OUT_DIR, "01_estimates_all.csv"),       row.names = FALSE)
write.csv(cmp,          file.path(OUT_DIR, "02_comparison_with_mimic.csv"), row.names = FALSE)
write.csv(interactions, file.path(OUT_DIR, "03_interactions.csv"),        row.names = FALSE)
if (!is.null(gradient)) write.csv(gradient, file.path(OUT_DIR, "04_gradient_within_cancer.csv"), row.names = FALSE)
write.csv(dist,         file.path(OUT_DIR, "05_trajectory_distribution.csv"), row.names = FALSE)
write.csv(dist_p,       file.path(OUT_DIR, "05b_trajectory_distribution_test.csv"), row.names = FALSE)
write.csv(balance,      file.path(OUT_DIR, "06_balance.csv"),             row.names = FALSE)

cat("\nDone. Aggregate results were written to", OUT_DIR, "\n")
print(cmp[, c("analysis", "rr", "rr_lo", "rr_hi", "rd", "rd_lo", "rd_hi")])
print(interactions)
