# =====================================================================
# 09_control_group_comparison.R  （第 2 版）
# 目的：加入"非癌症 AKI 对照组"，回答：
#   (1) AKI 轨迹分布在癌症和非癌症患者中是否不同？
#   (2) 轨迹与死亡的关系是否因癌症而改变？（交互检验）
#   (3) "出院后死亡负担"是不是癌症特有的？
#
# 运行前要有（都在 data_dir 里，csv 或 xlsx 均可）：
#   daily_stage.rds, traj_rule_based.rds   06 生成的（癌症组）
#   baseline_ctrl, outcomes_ctrl           04b、04e 把 cohort 换成 cohort_ctrl 后导出
#   daily_stage_ctrl                       04g 导出（对照组）
#   daily_stage_cancer_check（可选）       04g 换成 cohort 后导出，用来核对 SQL
# R 会话里最好还留着癌症组的 baseline、outc；没有的话脚本会弹窗让你选文件。
#
# 首次运行先装包（只需一次）：
# install.packages(c("dplyr","ggplot2","MatchIt","sandwich","lmtest","readr","readxl"))
#
# 【Excel 读取版】癌症组和对照组的文件无论扩展名是 .csv、.xlsx 还是被 Excel 改过，
#   脚本都会先识别真实格式再读取（read_any 用 readxl::excel_format 判断）。
#   可以在下面“文件路径”处直接填 5 个路径，填了就不会弹窗；留空（""）则弹窗选择。
# Data are read from data_dir (see the README).
# Data stay on the local computer; patient-level files must not be uploaded.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
  library(MatchIt); library(sandwich); library(lmtest)
})

if (!exists("data_dir")) data_dir <- getwd()   # folder with the exported CSV/RDS files
set.seed(20261019)

# ---- 文件路径（可选）：填了就直接读，留 "" 则弹窗让你选 ----
path_baseline_c <- ""    # 癌症组 baseline
path_outc_c     <- ""    # 癌症组 outcomes
path_baseline_k <- ""    # 对照组 baseline_ctrl
path_outc_k     <- ""    # 对照组 outcomes_ctrl
path_daily_k    <- ""    # 对照组 daily_stage_ctrl

# ---------------------------------------------------------------------
# 0. 读文件的小工具
# ---------------------------------------------------------------------
read_any <- function(path) {
  # 按文件内容而不是扩展名判断：Excel 另存的“.csv”其实可能是 xlsx
  fmt <- tryCatch(readxl::excel_format(path), error = function(e) NA_character_)
  if (!is.na(fmt)) {
    message("  读取 Excel 格式(", fmt, ")：", basename(path))
    d <- as.data.frame(readxl::read_excel(path))
  } else {
    message("  读取文本 csv：", basename(path))
    d <- as.data.frame(readr::read_csv(path, show_col_types = FALSE))
  }
  d
}
get_file <- function(path, stem, what) {   # 优先用填写的路径，其次在 data_dir 找 stem.csv/.xlsx，最后弹窗
  if (nzchar(path) && file.exists(path)) return(path)
  hit <- file.path(data_dir, paste0(stem, c(".csv", ".xlsx")))
  hit <- hit[file.exists(hit)]
  if (length(hit)) return(hit[1])

  stop("File not found: ", what, ". Set data_dir or the path_* variables.")
}
pick <- function(path, what) {            # 有路径就用路径，没有就弹窗
  if (nzchar(path) && file.exists(path)) return(path)

  stop("File not found: ", what, ". Set the path_* variables.")
}
find_file <- function(stem, required = TRUE) {
  cand <- file.path(data_dir, paste0(stem, c(".csv", ".xlsx")))
  hit  <- cand[file.exists(cand)]
  if (length(hit)) return(hit[1])
  if (!required) return(NA_character_)

  stop("File not found: ", stem, " in ", data_dir)
}
fix_outc <- function(o) {
  for (v in c("any_rrt", "death_28d", "death_90d")) o[[v]] <- as.logical(o[[v]])
  o
}

# ---------------------------------------------------------------------
# 1. 06 脚本的分型规则，原样封装成函数（癌症组和对照组必须用同一个）
#    输入：每人每天最高分期 daily(stay_id, day, max_stage)、baseline、outc
# ---------------------------------------------------------------------
safe_max <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }
to_time  <- function(x) if (inherits(x, "POSIXt")) x else as.POSIXct(as.character(x), tz = "UTC")

make_traj <- function(daily, baseline, outc) {
  traj <- daily %>%
    group_by(stay_id) %>%
    summarise(
      day0_stage  = safe_max(max_stage[day == 0]),
      obs_1_3     = sum(day %in% 1:3 & !is.na(max_stage)),
      persist_1_3 = sum(day %in% 1:3 & !is.na(max_stage) & max_stage >= 2),
      peak_later  = safe_max(max_stage[day >= 1]),
      late_aki    = any(day %in% 4:6 & !is.na(max_stage) & max_stage >= 2),
      .groups = "drop")

  baseline %>%
    select(stay_id, heme, solid, admission_age, sofa_day1) %>%
    left_join(outc %>% select(stay_id, t0, deathtime, hospital_expire_flag,
                              any_rrt, death_28d, death_90d), by = "stay_id") %>%
    left_join(traj, by = "stay_id") %>%
    mutate(
      obs_1_3     = ifelse(is.na(obs_1_3), 0, obs_1_3),
      persist_1_3 = ifelse(is.na(persist_1_3), 0, persist_1_3),
      late_aki    = ifelse(is.na(late_aki), FALSE, late_aki),
      death_h     = as.numeric(difftime(to_time(deathtime), to_time(t0), units = "hours"))
    ) %>%
    mutate(
      pheno = case_when(
        !is.na(death_h) & death_h <= 72                                   ~ "1_early_death_72h",
        obs_1_3 == 0                                                      ~ "2_no_data_day1to3",
        persist_1_3 == 0 & !late_aki                                      ~ "3_transient",
        persist_1_3 == 0 &  late_aki                                      ~ "4_recurrent_late",
        !is.na(peak_later) & peak_later >= 3 & (is.na(day0_stage) | day0_stage < 3) ~ "6_worsening_to_stage3",
        !is.na(day0_stage) & day0_stage >= 3                              ~ "7_severe_from_start",
        TRUE                                                              ~ "5_persistent_stage2"),
      tumor = case_when(heme == 1 ~ "heme", solid == 1 ~ "solid_only", TRUE ~ "none"))
}

# ---------------------------------------------------------------------
# 2. 癌症组：读入，并用函数重新分型，和 06 的结果核对
# ---------------------------------------------------------------------
# 注意：每次都重新读文件，不再沿用会话里可能残留的旧 baseline / outc
baseline <- read_any(pick(path_baseline_c, "癌症组 baseline 文件"))
outc     <- read_any(pick(path_outc_c,     "癌症组 outcomes 文件"))
cat("癌症组 baseline 行数:", nrow(baseline), " outcomes 行数:", nrow(outc), "（都应为 4213）\n")
need_b <- c("stay_id","gender","sepsis3","heme","solid","admission_age","sofa_day1")
need_o <- c("stay_id","t0","deathtime","dod","hospital_expire_flag","any_rrt","death_28d","death_90d")
if (nrow(baseline) != 4213 || nrow(outc) != 4213) stop("癌症组行数不是 4213，文件可能选错或读错。")
if (!all(need_b %in% names(baseline))) stop("癌症组 baseline 缺列：", paste(setdiff(need_b, names(baseline)), collapse = ", "))
if (!all(need_o %in% names(outc)))     stop("癌症组 outcomes 缺列：", paste(setdiff(need_o, names(outc)), collapse = ", "))
baseline_c <- as.data.frame(baseline)
outc_c     <- fix_outc(as.data.frame(outc))
daily_c    <- readRDS(file.path(data_dir, "daily_stage.rds"))

dat_c <- make_traj(daily_c, baseline_c, outc_c)

expected <- c("1_early_death_72h" = 388, "2_no_data_day1to3" = 1167, "3_transient" = 865,
              "4_recurrent_late" = 50, "5_persistent_stage2" = 835,
              "6_worsening_to_stage3" = 334, "7_severe_from_start" = 574)
got <- as.numeric(table(factor(dat_c$pheno, levels = names(expected))))
cat("\n===== 2. 核对：封装后的规则能否复现 06 的癌症组结果 =====\n")
print(data.frame(pheno = names(expected), expected = expected, got = got, row.names = NULL))
if (!all(got == expected)) stop("Reproduction check failed: the trajectory counts do not match the expected values.")
cat("一致，规则封装正确。\n")

# ---------------------------------------------------------------------
# 3. 对照组：读入、检查、分型
# ---------------------------------------------------------------------
baseline_k <- read_any(get_file(path_baseline_k, "baseline_ctrl", "对照组 baseline_ctrl"))
outc_k     <- fix_outc(read_any(get_file(path_outc_k, "outcomes_ctrl", "对照组 outcomes_ctrl")))
daily_k    <- read_any(get_file(path_daily_k, "daily_stage_ctrl", "对照组 daily_stage_ctrl"))
if ("day_n" %in% names(daily_k)) daily_k <- rename(daily_k, day = day_n)
daily_k <- daily_k %>% mutate(stay_id = stay_id, day = as.numeric(day), max_stage = as.numeric(max_stage))

cat("\n===== 3. 对照组数据检查 =====\n")
cat("baseline_ctrl 行数:", nrow(baseline_k), " 不重复 stay_id:", n_distinct(baseline_k$stay_id), "\n")
cat("outcomes_ctrl 行数:", nrow(outc_k), " 不重复 stay_id:", n_distinct(outc_k$stay_id), "\n")
cat("daily_stage_ctrl 行数:", nrow(daily_k), " 覆盖人数:", n_distinct(daily_k$stay_id), "\n")
if (anyDuplicated(baseline_k$stay_id) || anyDuplicated(outc_k$stay_id)) stop("对照组 stay_id 有重复")
if (length(intersect(baseline_k$stay_id, baseline_c$stay_id)) > 0) stop("对照组和癌症组有重叠的 stay_id")
cat("对照组 heme/solid 合计（应为 0）:", sum(baseline_k$heme, na.rm = TRUE) + sum(baseline_k$solid, na.rm = TRUE), "\n")
cat("日分期范围 day:", paste(range(daily_k$day), collapse = "~"), " max_stage:",
    paste(range(daily_k$max_stage, na.rm = TRUE), collapse = "~"), "\n")

# 可选：核对 04g 的 SQL 和 R 里算出来的是否逐行一致
chk_path <- find_file("daily_stage_cancer_check", required = FALSE)
if (!is.na(chk_path)) {
  chk <- read_any(chk_path)
  if ("day_n" %in% names(chk)) chk <- rename(chk, day = day_n)
  a <- daily_c %>% transmute(stay_id, day = as.numeric(day), r_stage = as.numeric(max_stage))
  b <- chk     %>% transmute(stay_id, day = as.numeric(day), sql_stage = as.numeric(max_stage))
  j <- full_join(a, b, by = c("stay_id", "day"))
  n_diff <- sum(is.na(j$r_stage) != is.na(j$sql_stage) |
                (!is.na(j$r_stage) & !is.na(j$sql_stage) & j$r_stage != j$sql_stage))
  cat("\nSQL 与 R 的每日分期逐行比对：R 行数", nrow(a), " SQL 行数", nrow(b),
      " 不一致行数", n_diff, "（应为 0）\n")
  if (n_diff > 0) stop("The daily stage from the SQL export and from the R derivation disagree; the export cannot be used for the control group.")
} else {
  cat("\n（没找到 daily_stage_cancer_check，跳过 SQL 核对。强烈建议做一次。）\n")
}

dat_k <- make_traj(daily_k, as.data.frame(baseline_k), outc_k)

cat("\n对照组轨迹人数:\n"); print(table(dat_k$pheno))

# ---------------------------------------------------------------------
# 4. 合并协变量，组成分析数据
# ---------------------------------------------------------------------
lev <- names(expected)
num01 <- function(x) as.numeric(as.logical(x))
z0    <- function(x) ifelse(is.na(x), 0, x)      # 合并症/脓毒症缺失 -> 视为"没有"（假设）

cov_cols <- c("gender", "sepsis3", "congestive_heart_failure", "renal_disease",
              "chronic_pulmonary_disease", "diabetes_without_cc", "diabetes_with_cc",
              "mild_liver_disease", "severe_liver_disease")
prep <- function(dat, base, cancer) {
  miss <- setdiff(cov_cols, names(base))
  if (length(miss)) stop("baseline 缺少这些列：", paste(miss, collapse = ", "))
  d <- dat %>% left_join(base %>% select(stay_id, all_of(cov_cols)), by = "stay_id")
  d %>% transmute(
    stay_id, cancer = cancer, group = ifelse(cancer == 1, "Cancer", "Non-cancer"),
    pheno = factor(pheno, levels = lev),
    age = admission_age, male = as.numeric(gender == "M"), sofa = sofa_day1,
    sepsis = z0(num01(sepsis3)),
    chf  = z0(num01(congestive_heart_failure)),
    ckd  = z0(num01(renal_disease)),
    copd = z0(num01(chronic_pulmonary_disease)),
    dm   = pmax(z0(num01(diabetes_without_cc)), z0(num01(diabetes_with_cc))),
    liver = pmax(z0(num01(mild_liver_disease)), z0(num01(severe_liver_disease))),
    hosp_death = num01(hospital_expire_flag), d28 = num01(death_28d),
    d90 = num01(death_90d), rrt = num01(any_rrt))
}

cat("\n性别编码（脚本假设 M/F）：\n")
print(table(baseline_c$gender, useNA = "ifany")); print(table(baseline_k$gender, useNA = "ifany"))

all_raw <- bind_rows(prep(dat_c, baseline_c, 1), prep(dat_k, as.data.frame(baseline_k), 0))
covs <- c("age", "male", "sofa", "sepsis", "chf", "ckd", "copd", "dm", "liver")
all_cc <- all_raw %>% filter(if_all(all_of(covs), ~ !is.na(.x)))

cat("\n===== 4. 样本量 =====\n")
cat("癌症组:", sum(all_raw$cancer == 1), " 对照组:", sum(all_raw$cancer == 0), "\n")
cat("协变量完整后 -> 癌症:", sum(all_cc$cancer == 1), " 对照:", sum(all_cc$cancer == 0),
    " (丢失", nrow(all_raw) - nrow(all_cc), "人)\n")

# ---------------------------------------------------------------------
# 5. 倾向评分匹配（1:2，最近邻，不放回，卡钳 0.2 SD）
#    只匹配基线变量，不匹配轨迹（轨迹是结局的一部分）
# ---------------------------------------------------------------------
m  <- matchit(reformulate(covs, response = "cancer"), data = all_cc,
              method = "nearest", distance = "glm", ratio = 2,
              caliper = 0.2, replace = FALSE)
md <- match.data(m)

cat("\n===== 5. 匹配结果 =====\n")
cat("匹配上的癌症患者:", sum(md$cancer == 1), "/", sum(all_cc$cancer == 1),
    " 匹配上的对照:", sum(md$cancer == 0), "\n")
smd <- function(x, g) {
  s <- sqrt((var(x[g == 1]) + var(x[g == 0])) / 2)
  (mean(x[g == 1]) - mean(x[g == 0])) / s
}
bal <- data.frame(variable = covs,
                  SMD_before = sapply(covs, function(v) smd(all_cc[[v]], all_cc$cancer)),
                  SMD_after  = sapply(covs, function(v) smd(md[[v]], md$cancer)),
                  row.names = NULL)
print(transform(bal, SMD_before = round(SMD_before, 3), SMD_after = round(SMD_after, 3)))
cat("(匹配后 |SMD| 一般要 < 0.1)\n")
write.csv(bal, file.path(data_dir, "table_balance_smd.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 6. 轨迹分布：癌症 vs 匹配对照
# ---------------------------------------------------------------------
cat("\n===== 6. 轨迹分布（匹配后，列百分比） =====\n")
tab6 <- table(md$pheno, md$group)
print(round(100 * prop.table(tab6, margin = 2), 1))
cat("卡方检验 p =", signif(chisq.test(tab6)$p.value, 3), "（未考虑配对，仅作描述）\n")

pd <- md %>% count(group, pheno) %>% group_by(group) %>% mutate(pct = 100 * n / sum(n))
write.csv(pd, file.path(data_dir, "table_pheno_by_group.csv"), row.names = FALSE)
p3 <- ggplot(pd, aes(x = group, y = pct, fill = pheno)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(pct >= 4, sprintf("%.1f", pct), "")),
            position = position_stack(vjust = 0.5), size = 3) +
  labs(x = NULL, y = "Patients (%)", fill = "Trajectory",
       title = "AKI trajectory distribution: cancer vs matched non-cancer") +
  theme_minimal(base_size = 11)
ggsave(file.path(data_dir, "fig3_pheno_cancer_vs_control.png"), p3, width = 7, height = 5, dpi = 300)

# ---------------------------------------------------------------------
# 7. 各轨迹 × 分组的结局（匹配后）；late_death = 住院存活者的 90 天死亡
# ---------------------------------------------------------------------
cat("\n===== 7. 各轨迹 × 分组的结局（%） =====\n")
out_tab <- md %>%
  group_by(pheno, group) %>%
  summarise(n = n(),
            hosp_death = round(100 * mean(hosp_death, na.rm = TRUE), 1),
            death_90d  = round(100 * mean(d90, na.rm = TRUE), 1),
            late_death = round(100 * mean(d90[hosp_death == 0], na.rm = TRUE), 1),
            rrt        = round(100 * mean(rrt, na.rm = TRUE), 1),
            .groups = "drop")
print(as.data.frame(out_tab))
write.csv(out_tab, file.path(data_dir, "table_outcomes_by_pheno_group.csv"), row.names = FALSE)

# ---------------------------------------------------------------------
# 8. 界标回归（排除 72h 内死亡者，参照 = 短暂型）
#    标准误按配对集聚类；交互检验 = Wald 检验
# ---------------------------------------------------------------------
prep_lm <- function(d) {
  d <- d %>% filter(pheno != "1_early_death_72h")
  d$pheno <- relevel(droplevels(d$pheno), ref = "3_transient")
  d
}
cl_of <- function(d) if ("subclass" %in% names(d)) d$subclass else seq_len(nrow(d))
est <- function(fit, term, cl) {
  b  <- coef(fit)[term]
  se <- sqrt(diag(vcovCL(fit, cluster = cl)))[term]
  data.frame(term = term, OR = exp(b), lo = exp(b - 1.96 * se), hi = exp(b + 1.96 * se),
             p = 2 * pnorm(-abs(b / se)), row.names = NULL)
}
fmt <- function(df) { num <- sapply(df, is.numeric); df[num] <- lapply(df[num], signif, 3); df }

run_models <- function(d, y, label) {
  d  <- d[!is.na(d[[y]]), ]
  cl <- cl_of(d)
  cv <- paste(covs, collapse = " + ")
  m_add <- glm(as.formula(paste(y, "~ cancer + pheno +", cv)), data = d, family = binomial)
  m_int <- glm(as.formula(paste(y, "~ cancer * pheno +", cv)), data = d, family = binomial)
  wt <- waldtest(m_add, m_int, vcov = function(x) vcovCL(x, cluster = cl), test = "Chisq")
  cat("\n--- ", label, "（n =", nrow(d), "）---\n")
  cat("癌症 vs 非癌症 校正 OR（控制轨迹后）:\n"); print(fmt(est(m_add, "cancer", cl)), row.names = FALSE)
  cat("癌症 × 轨迹 交互作用 Wald 检验 p =", signif(wt$`Pr(>Chisq)`[2], 3), "\n")

  cat("按轨迹分层：癌症 vs 非癌症 校正 OR\n")
  rows <- list()
  for (p in levels(d$pheno)) {
    s <- d[d$pheno == p, ]
    if (length(unique(s$cancer)) < 2 || sum(s[[y]]) < 10 || sum(1 - s[[y]]) < 10) next
    f <- glm(as.formula(paste(y, "~ cancer +", cv)), data = s, family = binomial)
    r <- est(f, "cancer", cl_of(s)); r$pheno <- p; r$n <- nrow(s); rows[[p]] <- r
  }
  if (length(rows)) print(fmt(do.call(rbind, rows)[, c("pheno", "n", "OR", "lo", "hi", "p")]), row.names = FALSE)

  cat("轨迹梯度（参照=短暂型），分组内校正 OR\n")
  for (g in c(1, 0)) {
    s <- d[d$cancer == g, ]
    f <- glm(as.formula(paste(y, "~ pheno +", cv)), data = s, family = binomial)
    tt <- grep("^pheno", names(coef(f)), value = TRUE)
    cat(ifelse(g == 1, "[癌症组]\n", "[非癌症组]\n"))
    print(fmt(do.call(rbind, lapply(tt, function(t) est(f, t, cl_of(s))))), row.names = FALSE)
  }
  invisible(NULL)
}

cat("\n===== 8. 界标回归（匹配样本，排除 72h 内死亡） =====\n")
dl <- prep_lm(md)
cat("界标样本: 癌症", sum(dl$cancer == 1), " 对照", sum(dl$cancer == 0), "\n")
run_models(dl, "hosp_death", "结局：住院死亡")
run_models(dl, "d90",        "结局：90 天死亡")

# ---------------------------------------------------------------------
# 9. 出院后死亡负担：只看住院存活者的 90 天死亡
# ---------------------------------------------------------------------
cat("\n===== 9. 住院存活者的 90 天死亡 =====\n")
surv <- dl %>% filter(hosp_death == 0)
print(as.data.frame(surv %>% group_by(group) %>%
        summarise(n = n(), death_90d_pct = round(100 * mean(d90, na.rm = TRUE), 1), .groups = "drop")))
run_models(surv, "d90", "结局：住院存活者的 90 天死亡")

# ---------------------------------------------------------------------
# 10. 敏感性：不匹配，全部癌症 + 全部对照，直接校正
# ---------------------------------------------------------------------
cat("\n===== 10. 敏感性分析：未匹配全样本（校正） =====\n")
dl_all <- prep_lm(all_cc)
cat("样本: 癌症", sum(dl_all$cancer == 1), " 对照", sum(dl_all$cancer == 0), "\n")
run_models(dl_all, "hosp_death", "未匹配：住院死亡")
run_models(dl_all, "d90",        "未匹配：90 天死亡")

cat("\n已保存到", data_dir, "：table_balance_smd.csv, table_pheno_by_group.csv,",
    "table_outcomes_by_pheno_group.csv, fig3_pheno_cancer_vs_control.png\n")
cat("Done.\n")


# ---------------------------------------------------------------------
# 11. 给后面的 07 / 08 / 10 脚本准备对象（它们按这些名字找数据）
# ---------------------------------------------------------------------
daily <- daily_c; dat <- dat_c; baseline <- baseline_c; outc <- outc_c
cat("已准备好 daily / dat / baseline / outc，可以继续运行 07、08、10。\n")
