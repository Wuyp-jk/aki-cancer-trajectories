make_traj6 <-
function(daily, outc) {
  tr <- daily |>
    group_by(stay_id) |>
    summarise(
      day0 = suppressWarnings(max(max_stage[day == 0], na.rm = TRUE)),
      n12  = sum(day %in% 1:2 & !is.na(max_stage)),
      ge2  = sum(day %in% 1:2 & !is.na(max_stage) & max_stage >= 2),
      s3   = sum(day %in% 1:2 & !is.na(max_stage) & max_stage >= 3),
      .groups = "drop") |>
    mutate(day0 = ifelse(is.infinite(day0), NA, day0))

  outc |>
    select(stay_id, t0, deathtime, dod, hospital_expire_flag) |>
    left_join(tr, by = "stay_id") |>
    mutate(
      n12 = coalesce(n12, 0L), ge2 = coalesce(ge2, 0L), s3 = coalesce(s3, 0L),
      death_h = as.numeric(difftime(to_time(deathtime), to_time(t0), units = "hours")),
      gap_d   = as.numeric(as.Date(substr(as.character(dod), 1, 10)) -
                           as.Date(substr(as.character(t0),  1, 10))),
      early_death = (!is.na(death_h) & death_h <= 72) |
                    (is.na(deathtime) & !is.na(gap_d) & gap_d <= 2),
      pheno6 = case_when(
        early_death                                  ~ "1_early_death",
        n12 == 0                                     ~ "2_early_icu_exit",
        ge2 == 0                                     ~ "3_transient",
        s3 > 0 & (is.na(day0) | day0 < 3)            ~ "5_worsening",
        !is.na(day0) & day0 >= 3                     ~ "6_severe_from_start",
        TRUE                                         ~ "4_persistent_stage2")) |>
    select(stay_id, pheno6, early_death, death_h, gap_d)
}
to_time <-
function(x) {
  if (inherits(x, "POSIXt")) return(x)
  x <- as.character(x)
  out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  i <- is.na(out) & !is.na(x) & nzchar(x)
  if (any(i)) out[i] <- as.POSIXct(x[i], tz = "UTC", format = "%Y-%m-%d")
  out
}
or_row <-
function(dat, outcome, extra = "") {
  f   <- glm(as.formula(paste(outcome, "~ cancer +", covs, extra)),
             data = dat, family = binomial)
  vc  <- vcovCL(f, cluster = dat$subclass)
  est <- coef(f)["cancer"]; se <- sqrt(vc["cancer", "cancer"])
  c(n = nrow(dat), OR = unname(exp(est)),
    lo = unname(exp(est - 1.96 * se)), hi = unname(exp(est + 1.96 * se)))
}
kid_fun <-
function(outc, base) {
  outc |>
    select(stay_id, dischtime, last_creat, last_creat_time, any_rrt, rrt_h, hospital_expire_flag) |>
    left_join(base |> select(stay_id, creat_base = creat_baseline_7d), by = "stay_id") |>
    mutate(
      rrt_any   = as.logical(any_rrt) %in% TRUE,
      gap_dc    = as.numeric(difftime(to_time(dischtime), to_time(last_creat_time), units = "hours")),
      rrt_early = rrt_any & rrt_h <= 72,
      survivor  = hospital_expire_flag == 0,
      evaluable = survivor & !is.na(creat_base) & creat_base >= 0.2 &
                  !is.na(last_creat) & !is.na(gap_dc) & gap_dc <= 72,
      y_rrt_any  = as.numeric(rrt_any),
      y_rrt_late = ifelse(rrt_early, NA, as.numeric(rrt_any & rrt_h > 72)),
      y_nonrec   = ifelse(evaluable, as.numeric(rrt_any | last_creat > 1.5 * creat_base), NA),
      y_make     = as.numeric(hospital_expire_flag == 1 | rrt_any |
                              (!is.na(last_creat) & !is.na(creat_base) & last_creat >= 2 * creat_base))) |>
    select(stay_id, starts_with("y_"), evaluable)
}
