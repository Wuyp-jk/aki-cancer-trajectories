# ============================================================
# 07  匹配队列的轨迹分析（更新版，严格72小时口径）
# 取代原来的 07_followup_checks.r（那份把早死类也放进模型，出现完全分离，OR 达 1e8~1e9，不能使用）。
#
# 前提：R 会话里已有 09 脚本留下的对象：
#       md（匹配后数据，含 subclass）、all_cc（匹配前、协变量完整的全体）、
#       daily_c, daily_k（每日最高分期）、dat_c, dat_k（含 pheno, day0_stage）
#       并已定义 covs、safe_max。没有时本脚本会自己补上。
#
# 分类：只用第0~2天（0~71小时）；早死 = t0后72小时内死亡（含只有死亡日期的院外死亡，见 early_ids）；
#       复发型并入短暂型；分析从72小时界标开始，早死类不进入死亡率模型。
# Output: prints summary numbers only; also saves 3 summary csv files (no patient-level data).
#       同时在工作目录保存 3 个汇总 csv（只有汇总数字）。
# ============================================================

library(dplyr)
library(sandwich)
library(lmtest)

need <- c("md", "all_cc", "daily_c", "daily_k", "dat_c", "dat_k", "outc_c", "outc_k")
miss <- need[!sapply(need, exists)]
if (length(miss) > 0) stop("缺少对象：", paste(miss, collapse = ", "), "。请先运行 09 脚本里的读数、分型和匹配部分。")

if (!exists("safe_max")) safe_max <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }
covs <- c("age", "male", "sofa", "sepsis", "chf", "ckd", "copd", "dm", "liver")
cov_str <- paste(covs, collapse = " + ")

# ---------- 1. 严格72小时分类（只用第0、1、2天） ----------
# 早死（72小时内死亡）的定义：
#   院内死亡有精确时间(death_h) -> death_h <= 72；
#   没有院内死亡时间的院外死亡，只有死亡日期 dod -> 日期差 <= 2 天（最晚也在 72 小时内，可以确定）；
#   日期差 = 3 天的死亡发生在 48~96 小时之间，无法确定，按界标之后处理（保守）。
to_time <- function(x) if (inherits(x, "POSIXt")) x else as.POSIXct(as.character(x), tz = "UTC")
if (!"dod" %in% names(outc_k)) stop("The control-group outcomes table has no dod column; check names(outc_k).")
early_ids <- function(dat, outc) {
  d <- dat %>% select(stay_id, death_h) %>%
    left_join(outc %>% select(stay_id, dod, t0), by = "stay_id") %>%
    mutate(dod_days = as.numeric(as.Date(to_time(dod)) - as.Date(to_time(t0))),
           early72  = (!is.na(death_h) & death_h <= 72) |
                      (is.na(death_h) & !is.na(dod_days) & dod_days <= 2))
  d$stay_id[d$early72]
}
early_c <- early_ids(dat_c, outc_c); early_k <- early_ids(dat_k, outc_k)
cat("\n早死人数：癌症组", length(early_c), "（其中新增的院外死亡", length(early_c) - sum(dat_c$pheno == "1_early_death_72h"),
    "人）；对照组", length(early_k), "（新增", length(early_k) - sum(dat_k$pheno == "1_early_death_72h"), "人）\n")

re72_strict <- function(daily, dat, early) {
  w <- daily %>% group_by(stay_id) %>%
    summarise(
      peak_1_2    = safe_max(max_stage[day %in% 1:2]),
      persist_1_2 = sum(day %in% 1:2 & !is.na(max_stage) & max_stage >= 2),
      obs_1_2     = sum(day %in% 1:2 & !is.na(max_stage)),
      .groups = "drop")
  dat %>% left_join(w, by = "stay_id") %>%
    mutate(
      obs_1_2     = ifelse(is.na(obs_1_2), 0, obs_1_2),
      persist_1_2 = ifelse(is.na(persist_1_2), 0, persist_1_2),
      pheno72 = case_when(
        stay_id %in% early                                                       ~ "1_early_death_72h",
        obs_1_2 == 0                                                             ~ "2_no_data_day1to2",
        persist_1_2 == 0                                                         ~ "3_transient",
        !is.na(peak_1_2) & peak_1_2 >= 3 & (is.na(day0_stage) | day0_stage < 3)  ~ "6_worsening_to_stage3",
        !is.na(day0_stage) & day0_stage >= 3                                     ~ "7_severe_from_start",
        TRUE                                                                     ~ "5_persistent_stage2")) %>%
    select(stay_id, pheno72)
}

lev72 <- c("1_early_death_72h", "2_no_data_day1to2", "3_transient",
           "5_persistent_stage2", "6_worsening_to_stage3", "7_severe_from_start")
p72 <- bind_rows(re72_strict(daily_c, dat_c, early_c), re72_strict(daily_k, dat_k, early_k))

md_unified <- md %>% select(-any_of("pheno72")) %>% left_join(p72, by = "stay_id") %>%
  mutate(pheno72 = factor(pheno72, levels = lev72))

cat("\n===== A. 轨迹分布（匹配队列，严格72小时版本） =====\n")
tab3 <- table(md_unified$pheno72, md_unified$cancer)
print(tab3)
print(round(100 * prop.table(tab3, margin = 2), 1))
ch <- chisq.test(tab3)
cat("卡方 =", round(ch$statistic, 1), " 自由度 =", ch$parameter, " P =", signif(ch$p.value, 3), "\n")
cat("各类单独的2x2卡方 P 值：\n")
for (i in seq_len(nrow(tab3))) {
  t2 <- rbind(tab3[i, ], colSums(tab3) - tab3[i, ])
  cat(sprintf("  %-24s P = %.3g\n", rownames(tab3)[i], suppressWarnings(chisq.test(t2, correct = FALSE)$p.value)))
}
write.csv(as.data.frame.matrix(tab3), "table3_counts_72h.csv")

md_surv <- md_unified %>% filter(pheno72 != "1_early_death_72h")
md_surv$pheno72 <- droplevels(md_surv$pheno72)
cat("\n界标后人数:", nrow(md_surv), " 癌症:", sum(md_surv$cancer == 1), " 非癌症:", sum(md_surv$cancer == 0), "\n")

# ---------- 2. 总体 OR（两种口径都给，请告诉我论文原来用哪一种） ----------
cl_or <- function(f, data, term = "cancer") {
  b <- coef(f)[term]; se <- sqrt(vcovCL(f, cluster = data$subclass)[term, term])
  round(c(OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se)), 2)
}
cat("\n===== B. 癌症对死亡的总体 OR（早死已剔除，聚类稳健标准误） =====\n")
for (y in c("hosp_death", "d90")) {
  f1 <- glm(as.formula(paste(y, "~ cancer +", cov_str)), data = md_surv, family = binomial)
  f2 <- glm(as.formula(paste(y, "~ cancer + pheno72 +", cov_str)), data = md_surv, family = binomial)
  cat(y, "| 只校正协变量:", cl_or(f1, md_surv), "\n")
  cat(y, "| 再校正轨迹    :", cl_or(f2, md_surv), "\n")
}

# ---------- 3. 各轨迹内癌症对死亡的 OR（表4 / 图4 / 图S1） ----------
cat("\n===== C. 各轨迹内：癌症 vs 非癌症 的 OR =====\n")
res_list <- list()
for (p in levels(md_surv$pheno72)) {
  sub <- md_surv[md_surv$pheno72 == p, ]
  fh <- glm(as.formula(paste("hosp_death ~ cancer +", cov_str)), data = sub, family = binomial)
  f9 <- glm(as.formula(paste("d90 ~ cancer +", cov_str)), data = sub, family = binomial)
  h <- cl_or(fh, sub); d <- cl_or(f9, sub)
  res_list[[p]] <- data.frame(Trajectory = p, N_Cancer = sum(sub$cancer == 1), N_Control = sum(sub$cancer == 0),
    Hosp_OR = h[1], Hosp_Lo = h[2], Hosp_Hi = h[3], D90_OR = d[1], D90_Lo = d[2], D90_Hi = d[3])
}
res_df <- do.call(rbind, res_list); rownames(res_df) <- NULL
print(res_df)
write.csv(res_df, "table4_or_by_trajectory_72h.csv", row.names = FALSE)

# ---------- 4. 交互作用 ----------
cat("\n===== D. 癌症 × 轨迹 交互作用 Wald P 值 =====\n")
for (y in c("hosp_death", "d90")) {
  fa <- glm(as.formula(paste(y, "~ cancer + pheno72 +", cov_str)), data = md_surv, family = binomial)
  fi <- glm(as.formula(paste(y, "~ cancer * pheno72 +", cov_str)), data = md_surv, family = binomial)
  pv <- waldtest(fa, fi, test = "Chisq", vcov = function(x) vcovCL(x, cluster = md_surv$subclass))$`Pr(>Chisq)`[2]
  cat(y, ": P =", signif(pv, 3), "\n")
}

# ---------- 5. 各类粗死亡率（图5和正文的绝对差） ----------
cat("\n===== E. 各类粗死亡率 (%) =====\n")
crude <- md_surv %>% group_by(pheno72, cancer) %>%
  summarise(n = n(), hosp = round(100 * mean(hosp_death), 1), d90 = round(100 * mean(d90), 1), .groups = "drop") %>%
  as.data.frame()
print(crude)
write.csv(crude, "crude_mortality_by_trajectory_72h.csv", row.names = FALSE)

# ---------- 6. 癌症组 / 非癌症组各自：各轨迹相对短暂型的 OR ----------
cat("\n===== F. 组内：各轨迹相对短暂型的 OR =====\n")
for (g in c(1, 0)) {
  s <- md_surv[md_surv$cancer == g, ]
  s$pheno72 <- relevel(droplevels(s$pheno72), "3_transient")
  for (y in c("d90", "hosp_death")) {
    f <- glm(as.formula(paste(y, "~ pheno72 +", cov_str)), data = s, family = binomial)
    b <- coef(f); se <- sqrt(diag(vcovCL(f, cluster = s$subclass))); k <- grep("pheno72", names(b))
    cat("\n", ifelse(g == 1, "癌症组", "非癌症组"), "|", y, "\n")
    print(round(cbind(OR = exp(b[k]), lo = exp(b[k] - 1.96 * se[k]), hi = exp(b[k] + 1.96 * se[k])), 2))
  }
}

# ---------- 7. 未匹配队列的敏感性分析 ----------
cat("\n===== G. 未匹配队列（协变量校正，无聚类） =====\n")
all_u <- all_cc %>% select(-any_of("pheno72")) %>% left_join(p72, by = "stay_id") %>%
  filter(pheno72 != "1_early_death_72h") %>%
  mutate(pheno72 = droplevels(factor(pheno72, levels = lev72)))
cat("未匹配 n =", nrow(all_u), " 癌症:", sum(all_u$cancer == 1), " 非癌症:", sum(all_u$cancer == 0), "\n")
ci_cancer <- function(f) round(exp(cbind(OR = coef(f), confint.default(f)))["cancer", ], 2)
for (y in c("hosp_death", "d90")) {
  cat(y, "| 总体(协变量):", ci_cancer(glm(as.formula(paste(y, "~ cancer +", cov_str)), data = all_u, family = binomial)), "\n")
  cat(y, "| 总体(+轨迹)  :", ci_cancer(glm(as.formula(paste(y, "~ cancer + pheno72 +", cov_str)), data = all_u, family = binomial)), "\n")
}
cat("\n90天死亡，各轨迹内:\n")
for (p in levels(all_u$pheno72)) {
  s <- all_u[all_u$pheno72 == p, ]
  cat(sprintf("  %-24s", p), ci_cancer(glm(as.formula(paste("d90 ~ cancer +", cov_str)), data = s, family = binomial)), "\n")
}
fa <- glm(as.formula(paste("d90 ~ cancer + pheno72 +", cov_str)), data = all_u, family = binomial)
fi <- glm(as.formula(paste("d90 ~ cancer * pheno72 +", cov_str)), data = all_u, family = binomial)
cat("未匹配交互 P (90天) =", signif(anova(fa, fi, test = "Chisq")$`Pr(>Chi)`[2], 3), "\n")

cat("\n已保存汇总 csv 到：", getwd(), "\n  table3_counts_72h.csv\n  table4_or_by_trajectory_72h.csv\n  crude_mortality_by_trajectory_72h.csv\n")
