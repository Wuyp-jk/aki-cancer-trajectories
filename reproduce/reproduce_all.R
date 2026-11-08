# =============================================================================
# reproduce_all.R  (v8:主分析校正集改为 B = 9 个匹配协变量 + 入院类型/ICU 类型/入院科室/年份段)
# 基于 reproduce_all.R (v4:效应量为风险比 RR 与风险差 RD)。预期值已更新为 v8 稿子的数字(主校正集 B);
# 全部检查项应为 OK。结果表见 reproduce_output/ 及 v7_*.csv。
# Trajectories of acute kidney injury and 90-day mortality in critically ill
# patients with and without cancer (MIMIC-IV v3.1)
#
# 复现稿子(v8)中的 Table 1–7、补充图 S1、补充表 S1–S5 的数字、AUC,
# 以及 RRT / 安宁疗护 / 去肾脏成分 SOFA 的核查。
#
# 每个检查项都带状态标记:
#   [R]  2026-10-23 已在最终版数据上完全复现(Table 1、2、3 与 RRT 核查,沿用旧版)
#   [V4] v4 新增:预期值取自稿子定稿时的运行结果(RR / RD 的函数与本脚本相同)
#        v8 稿子(主校正集 B)在最终版数据上端到端运行(R 4.6.1):893 项检查全部 OK,0 项 DIFF;
#        第 13 节因缺 procedureevents 未运行(见下方运行前说明)
#   [V8] 第 15 节(无条件估计、逐步缩小样本、早期出 ICU 拆分、平衡表、E 值、code status 记录)的预期值取自 v8 稿子;
#        队列内 LRT / AUC 的预期值同样已更新为主校正集 B 的结果(折在合并匹配样本上分配)
#
# RR:修正泊松回归 + 稳健方差;RD:线性概率模型 + 稳健方差(单位:百分点)。
# 匹配样本按匹配集聚类;仅肿瘤患者、未匹配样本使用不聚类的稳健标准误。
#
# 运行前:
#   1. 在仍保留分析对象的会话里运行 00_export_inputs.R(生成 inputs_final.rds)
#   2. 设置 ADM_PATH(admissions.csv(.gz) 或自己导出的 hadm_id + discharge_location 两列 csv)
#      和 PE_PATH(procedureevents.csv(.gz));第 12、13 节需要。
#   3. 去肾脏成分 SOFA 的补充分析(第 12b 节)需要 SOFA_PATH:BigQuery 导出的 csv,
#      列为 stay_id, sofa, renal(来自 mimiciv_derived.first_day_sofa);没有就自动跳过。
# v7(主分析校正集 B):另需 EXTRA_PATH <- file.choose()  # extra_confounders.csv(sql/S14)
# 运行: source("reproduce_all.R", echo = FALSE, encoding = "UTF-8")
# 输出: reproduce_output/ 下的 checks.csv(逐项对照稿子)、各表 csv、session_info.txt
# =============================================================================

# ---- 0. 配置 -----------------------------------------------------------------
if (!exists("PROJ_DIR")) PROJ_DIR <- getwd()   # inputs_final.rds 所在文件夹;输出写入其下的 reproduce_output/
# 下面几个变量,如果在运行 source() 之前已经在控制台里设好,就沿用控制台里的值:
#   ADM_PATH <- file.choose()   # admissions.csv(.gz) 或 adm_hospice.csv
#   PE_PATH  <- file.choose()   # procedureevents.csv(.gz)
#   DI_PATH  <- file.choose()   # d_items.csv(.gz),可选;没有就设为 NA_character_
#   SOFA_PATH <- file.choose()  # sofa_renal.csv,可选
# 或者设 MIMIC_DIR(会在该文件夹内递归查找;不要设成主目录,否则查找会非常慢)
if (!exists("MIMIC_DIR")) MIMIC_DIR <- NA_character_
if (!exists("ADM_PATH"))  ADM_PATH  <- NA_character_   # 取 discharge_location
if (!exists("PE_PATH"))   PE_PATH   <- NA_character_   # procedureevents
if (!exists("DI_PATH"))   DI_PATH   <- NA_character_   # d_items(可选)
if (!exists("SOFA_PATH")) SOFA_PATH <- NA_character_   # stay_id, sofa, renal(可选)
OUT_DIR   <- file.path(PROJ_DIR, "reproduce_output")
dir.create(OUT_DIR, showWarnings = FALSE)
setwd(PROJ_DIR)

suppressPackageStartupMessages({
  library(dplyr); library(sandwich); library(lmtest); library(pROC)
})

# ---- 1. 读取输入 ---------------------------------------------------------------
inp <- readRDS(file.path(PROJ_DIR, "inputs_final.rds"))
daily_c <- inp$daily_c;  daily_k <- inp$daily_k
outc_c  <- inp$outc_c;   outc_k  <- inp$outc_k
baseline_c <- inp$baseline_c; baseline_k <- inp$baseline_k
md <- inp$md;  kdigo_temp <- inp$kdigo_temp
md$pheno <- NULL                                   # 旧版轨迹,不使用

cov_match <- c("age", "male", "sofa", "sepsis", "chf", "ckd", "copd", "dm", "liver")   # 匹配用的 9 个协变量(Table 2)

# ---- v7:主分析的校正集 B = 9 个匹配协变量 + 入院类型、首个 ICU 类型、入院科室、年份段 ----------
# 需要 EXTRA_PATH(sql/S14 导出的 extra_confounders.csv)。code status 只用于补充分析(第 15 节)。
if (!exists("EXTRA_PATH") || is.na(EXTRA_PATH) || !file.exists(EXTRA_PATH))
  stop("请先设置 EXTRA_PATH <- file.choose()   # extra_confounders.csv")
ex_raw <- read.csv(EXTRA_PATH, stringsAsFactors = FALSE)
stopifnot(all(c("stay_id", "cancer", "admission_type", "first_careunit", "anchor_year_group",
                "first_service", "code_status_raw") %in% names(ex_raw)), !anyDuplicated(ex_raw$stay_id))
ex_cov <- ex_raw |> mutate(
  adm_type = case_when(
    admission_type %in% c("ELECTIVE", "SURGICAL SAME DAY ADMISSION") ~ "elective/surgical",
    admission_type %in% c("EW EMER.", "DIRECT EMER.", "URGENT")      ~ "emergency/urgent",
    TRUE                                                              ~ "observation/other"),
  unit = case_when(
    grepl("MICU/SICU|Medical/Surgical", first_careunit)               ~ "MICU/SICU",
    grepl("Medical Intensive|\\(MICU\\)", first_careunit)             ~ "MICU",
    grepl("Cardiac|Coronary|CVICU|CCU", first_careunit)               ~ "Cardiac",
    grepl("Neuro", first_careunit)                                    ~ "Neuro",
    TRUE                                                              ~ "Surgical/trauma/other"),   # 含 SICU、TSICU 及 54 例其他单元
  period = ifelse(is.na(anchor_year_group) | anchor_year_group == "", "unknown", anchor_year_group),
  svc = ifelse(first_service %in% c("SURG", "CSURG", "TSURG", "VSURG", "NSURG", "ORTHO", "TRAUM",
                                    "PSURG", "ENT", "GU", "GYN"), "surgical", "medical/other"),
  cs_up = toupper(ifelse(is.na(code_status_raw), "", code_status_raw)),
  code = case_when(
    cs_up == ""                                                           ~ "not documented",
    grepl("COMFORT|CMO", cs_up)                                           ~ "comfort measures",
    grepl("DNR|DO NOT RESUS|DNI|DO NOT INTUB|NOT INDICATE|NO CPR", cs_up) ~ "DNR/DNI",
    grepl("FULL", cs_up)                                                  ~ "full code",
    TRUE                                                                  ~ "other")) |>
  select(stay_id, adm_type, unit, svc, period, code)
md <- md |> left_join(ex_cov, by = "stay_id")
stopifnot(!anyNA(md$adm_type), !anyNA(md$unit), !anyNA(md$period))
cov_extra <- c("adm_type", "unit", "svc", "period")
cov_names <- c(cov_match, cov_extra)                 # 模型校正集 B(主分析)
covs <- paste(cov_names, collapse = " + ")           # 公式字符串;向量形式用 cov_names
RES  <- list()                                     # 汇总各表,最后写 csv

# ---- 2. 工具函数 ---------------------------------------------------------------
# 时间转换:已是日期时间则原样返回;文字则按 "年-月-日 时:分:秒" 解析,
# 只有日期的字符串按日期解析,空字符串得 NA
to_time <- function(x) {
  if (inherits(x, "POSIXt")) return(x)
  x <- as.character(x)
  out <- as.POSIXct(x, tz = "UTC", format = "%Y-%m-%d %H:%M:%S")
  i <- is.na(out) & !is.na(x) & nzchar(x)
  if (any(i)) out[i] <- as.POSIXct(x[i], tz = "UTC", format = "%Y-%m-%d")
  out
}
to01 <- function(x) as.numeric(as.character(x) %in% c("1", "TRUE", "T"))
z0   <- function(x) ifelse(is.na(x), 0, x)
num  <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)
num01 <- to01

# 轨迹分类(稿子的 6 类规则:第 0–2 天,不使用晚期数据;
# 死亡日期不晚于 t0 后两个日历日、且无精确死亡时间者也算早期死亡)
make_traj6 <- function(daily, outc) {
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
    select(stay_id, t0, deathtime, dod) |>
    left_join(tr, by = "stay_id") |>
    mutate(
      n12 = coalesce(n12, 0L), ge2 = coalesce(ge2, 0L), s3 = coalesce(s3, 0L),
      death_h = as.numeric(difftime(to_time(deathtime), to_time(t0), units = "hours")),
      gap_d   = as.numeric(base::as.Date(substr(as.character(dod), 1, 10)) -
                           base::as.Date(substr(as.character(t0),  1, 10))),
      early_death = (!is.na(death_h) & death_h <= 72) |
                    (is.na(deathtime) & !is.na(gap_d) & gap_d <= 2),
      pheno6 = case_when(
        early_death                           ~ "1_early_death",
        n12 == 0                              ~ "2_early_icu_exit",
        ge2 == 0                              ~ "3_transient",
        s3 > 0 & (is.na(day0) | day0 < 3)     ~ "5_worsening",
        !is.na(day0) & day0 >= 3              ~ "6_severe_from_start",
        TRUE                                  ~ "4_persistent_stage2")) |>
    select(stay_id, pheno6, early_death, death_h, gap_d)
}

#<<helpers (RR / RD;与 aki_rr_rd.R 相同)
.prep <- function(data, outcome, terms, cluster) {
  vars <- unique(c(outcome, all.vars(reformulate(terms)), cluster))
  data[complete.cases(data[, vars, drop = FALSE]), , drop = FALSE]
}
.vcov <- function(m, cl) {
  if (is.null(cl)) vcovHC(m, type = "HC1") else vcovCL(m, cluster = cl, type = "HC1")
}
# RR:修正泊松回归
fit_rr <- function(data, outcome, terms, cluster = NULL) {
  d <- .prep(data, outcome, terms, cluster)
  f <- reformulate(terms, response = outcome)
  X <- model.matrix(f, d)
  m <- glm(f, data = d, family = poisson(link = "log"),
           start = c(log(mean(d[[outcome]])), rep(0, ncol(X) - 1)),
           control = glm.control(maxit = 100))
  cl <- if (is.null(cluster)) NULL else d[[cluster]]
  list(model = m, vcov = .vcov(m, cl), data = d, cl = cl, link = "log")
}
# RD:线性概率模型
fit_rd <- function(data, outcome, terms, cluster = NULL) {
  d <- .prep(data, outcome, terms, cluster)
  f <- reformulate(terms, response = outcome)
  m <- lm(f, data = d)
  cl <- if (is.null(cluster)) NULL else d[[cluster]]
  list(model = m, vcov = .vcov(m, cl), data = d, cl = cl, link = "identity")
}
tidy_fit <- function(fit, terms) {
  ct <- coeftest(fit$model, vcov. = fit$vcov)
  do.call(rbind, lapply(terms, function(tm) {
    b <- ct[tm, 1]; s <- ct[tm, 2]
    f <- if (fit$link == "log") exp else identity
    data.frame(term = tm, estimate = f(b), lo = f(b - 1.96 * s),
               hi = f(b + 1.96 * s), p = ct[tm, 4], row.names = NULL)
  }))
}
# 肿瘤 vs 非肿瘤:总体(含 by 变量)+ 各层;RR 与 RD 各一行
cancer_effect <- function(data, outcome, by, cluster = "subclass",
                          adj = cov_names, exposure = "cancer") {
  out <- list()
  for (sc in c("RR", "RD")) {
    fitter <- if (sc == "RR") fit_rr else fit_rd
    ov <- fitter(data, outcome, c(exposure, by, adj), cluster)
    res <- cbind(stratum = "All", scale = sc, tidy_fit(ov, exposure))
    for (lv in levels(factor(data[[by]]))) {
      d <- data[data[[by]] == lv, ]
      f <- fitter(d, outcome, c(exposure, adj), cluster)
      res <- rbind(res, cbind(stratum = lv, scale = sc, tidy_fit(f, exposure)))
    }
    out[[sc]] <- res
  }
  do.call(rbind, out)
}
# cancer × 分层变量的交互:乘法(RR)与加性(RD),稳健 Wald 检验
interaction_tests <- function(data, outcome, by, cluster = "subclass",
                              adj = cov_names, exposure = "cancer") {
  res <- list()
  for (sc in c("RR", "RD")) {
    fitter <- if (sc == "RR") fit_rr else fit_rd
    big <- fitter(data, outcome, c(paste0(exposure, "*", by), adj), cluster)
    sml <- fitter(data, outcome, c(exposure, by, adj), cluster)
    cl <- big$cl
    wt <- waldtest(sml$model, big$model, vcov = function(x) .vcov(x, cl), test = "Chisq")
    res[[sc]] <- data.frame(scale = sc, df = wt$Df[2], chisq = wt$Chisq[2],
                            p = wt$`Pr(>Chisq)`[2])
  }
  do.call(rbind, res)
}
# 仅肿瘤患者(独立观测):各组相对参照组
group_effects <- function(data, outcome, group, adj = cov_names, cluster = NULL) {
  data[[group]] <- factor(data[[group]])
  g_terms <- paste0(group, levels(data[[group]])[-1])
  res <- NULL
  for (sc in c("RR", "RD")) {
    fitter <- if (sc == "RR") fit_rr else fit_rd
    f <- fitter(data, outcome, c(group, adj), cluster)
    res <- rbind(res, cbind(scale = sc, tidy_fit(f, g_terms)))
  }
  res
}
crude_risk <- function(data, outcome, by, exposure = "cancer") {
  agg <- aggregate(data[[outcome]],
                   by = list(stratum = data[[by]], exposure = data[[exposure]]),
                   FUN = function(x) c(n = length(x), risk = mean(x)))
  agg <- do.call(data.frame, agg); names(agg)[3:4] <- c("n", "risk"); agg
}
#>>helpers

# 肾脏结局(口径按稿子的文字描述写成)
kid_fun <- function(outc, base) {
  outc |>
    select(stay_id, dischtime, last_creat, last_creat_time, any_rrt, rrt_h,
           hospital_expire_flag) |>
    left_join(base |> select(stay_id, creat_base = creat_baseline_7d), by = "stay_id") |>
    mutate(
      last_creat = as.numeric(last_creat), rrt_h = as.numeric(rrt_h),
      rrt_any   = as.logical(any_rrt) %in% TRUE,
      gap_dc    = as.numeric(difftime(to_time(dischtime), to_time(last_creat_time),
                                      units = "hours")),
      rrt_early = rrt_any & rrt_h <= 72,                 # "72 小时内"含恰好 72 小时
      survivor  = hospital_expire_flag == 0,
      evaluable = survivor & !is.na(creat_base) & creat_base >= 0.2 &
                  !is.na(last_creat) & !is.na(gap_dc) & gap_dc <= 72,
      y_rrt_any  = as.numeric(rrt_any),
      y_rrt_late = ifelse(rrt_early, NA, as.numeric(rrt_any & rrt_h > 72)),
      y_nonrec   = ifelse(evaluable, as.numeric(rrt_any | last_creat > 1.5 * creat_base), NA),
      y_make     = as.numeric(hospital_expire_flag == 1 | rrt_any |
                              (!is.na(last_creat) & !is.na(creat_base) &
                               last_creat >= 2 * creat_base))) |>
    select(stay_id, starts_with("y_"), evaluable)
}
wald_p <- function(f0, f1, cl) {
  w <- waldtest(f0, f1, vcov = function(x) vcovCL(x, cluster = cl), test = "Chisq")
  c(chisq = w$Chisq[2], df = w$Df[2], p = w[["Pr(>Chisq)"]][2])
}

# 读取 gz 压缩 csv 的指定列(只读需要的列,省内存;未压缩的 csv 也能读)
read_gz_cols <- function(path, cols) {
  hdr <- gsub('"', "", strsplit(readLines(gzfile(path), n = 1), ",")[[1]])
  cc  <- ifelse(hdr %in% cols, "character", "NULL"); names(cc) <- hdr
  read.csv(gzfile(path), colClasses = cc, na.strings = c("", "NA"))
}

# 逐项对照稿子
CHK <- NULL
chk <- function(item, observed, expected, tol = 0.0051, note = "") {
  obs <- if (length(observed) == 1) as.numeric(observed) else NA_real_
  ok  <- is.finite(obs) && abs(obs - expected) <= tol
  CHK <<- rbind(CHK, data.frame(item = item, expected = expected,
                                observed = round(obs, 4),
                                status = if (ok) "OK" else "DIFF", note = note,
                                stringsAsFactors = FALSE))
  invisible(ok)
}
# e = c(RR, 下限, 上限, RD, 下限, 上限);RD 以百分点计
chk_vec <- function(prefix, e, rr, rd, note = "") {
  chk(paste(prefix, "RR"),     rr$estimate,       e[1], 0.0051, note)
  chk(paste(prefix, "RR lo"),  rr$lo,             e[2], 0.0051, note)
  chk(paste(prefix, "RR hi"),  rr$hi,             e[3], 0.0051, note)
  chk(paste(prefix, "RD"),     100 * rd$estimate, e[4], 0.051,  note)
  chk(paste(prefix, "RD lo"),  100 * rd$lo,       e[5], 0.051,  note)
  chk(paste(prefix, "RD hi"),  100 * rd$hi,       e[6], 0.051,  note)
}
chk_eff <- function(prefix, res, exp, note = "") {          # cancer_effect() 的结果
  for (s in names(exp)) chk_vec(paste(prefix, s), exp[[s]],
    res[res$stratum == s & res$scale == "RR", ], res[res$stratum == s & res$scale == "RD", ], note)
}
chk_grp <- function(prefix, res, exp, note = "") {          # group_effects() 的结果(组变量 g5)
  for (s in names(exp)) chk_vec(paste(prefix, s), exp[[s]],
    res[res$scale == "RR" & res$term == paste0("g5", s), ],
    res[res$scale == "RD" & res$term == paste0("g5", s), ], note)
}
chk_terms <- function(prefix, rr, rd, exp, note = "") {     # tidy_fit() 的结果,按 term
  for (s in names(exp)) chk_vec(paste(prefix, s), exp[[s]],
    rr[rr$term == s, ], rd[rd$term == s, ], note)
}
chk_int <- function(prefix, tab, rr_chisq, rd_chisq) {      # interaction_tests() 的结果
  chk(paste(prefix, "乘法 chisq"), tab$chisq[tab$scale == "RR"], rr_chisq, tol = 0.01)
  chk(paste(prefix, "加性 chisq"), tab$chisq[tab$scale == "RD"], rd_chisq, tol = 0.01)
}

# ---- 3. 轨迹分类;Table 3 ---------------------------------------------------- [R]
traj_c <- make_traj6(daily_c, outc_c)
traj_k <- make_traj6(daily_k, outc_k)
traj_all <- bind_rows(traj_c, traj_k)

tc <- table(traj_c$pheno6)
for (i in 1:6) chk(paste("肿瘤队列", names(tc)[i]), as.numeric(tc[i]),
                   c(399, 1158, 989, 866, 235, 566)[i], tol = 0)

md2 <- md |> select(stay_id, cancer, subclass) |>
  left_join(traj_all |> select(stay_id, pheno6), by = "stay_id")
stopifnot(!anyNA(md2$pheno6))
t3 <- table(md2$cancer, md2$pheno6)
exp3 <- rbind(`0` = c(466, 2644, 2194, 1760, 385, 975),
              `1` = c(399, 1157,  989,  866, 235, 566))
for (g in c("0", "1")) for (i in 1:6)
  chk(paste("Table 3 cancer =", g, colnames(t3)[i]), t3[g, i], exp3[g, i], tol = 0)
RES$table3 <- as.data.frame.matrix(t3)

# ---- 4. Table 2:匹配后平衡(SMD)------------------------------------------- [R]
smd_fun <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
}
t2 <- data.frame(variable = cov_match,
                 smd = sapply(cov_match, function(v) smd_fun(num(md[[v]]), md$cancer)))
exp_smd <- c(age = .015, male = .004, sofa = .044, sepsis = .001, chf = .040,
             ckd = .038, copd = .033, dm = .026, liver = .038)
for (v in cov_match) chk(paste("Table 2 SMD", v), abs(t2$smd[t2$variable == v]),
                         exp_smd[[v]], tol = 0.0006)
RES$table2 <- t2

# ---- 5. 主分析样本;Table 4 与交互检验 ---------------------------------------- [V4]
four <- c("3_transient", "4_persistent_stage2", "5_worsening", "6_severe_from_start")
tlab <- c("Transient", "Persistent2", "Worsening", "Severe")
md3 <- md |>
  select(stay_id, cancer, subclass, all_of(cov_names), hosp_death, d90) |>
  left_join(traj_all |> select(stay_id, pheno6), by = "stay_id") |>
  filter(!pheno6 %in% c("1_early_death", "2_early_icu_exit"))
md3$traj <- factor(md3$pheno6, levels = four, labels = tlab)
chk("主样本 n", nrow(md3), 7970, tol = 0)
chk("主样本 肿瘤 n", sum(md3$cancer == 1), 2656, tol = 0)
chk("主样本 非肿瘤 n", sum(md3$cancer == 0), 5314, tol = 0)

# 预期值 c(RR, lo, hi, RD, lo, hi);RD 为百分点
t4_90 <- cancer_effect(md3, "d90", "traj")
t4_h  <- cancer_effect(md3, "hosp_death", "traj")
chk_eff("Table 4 d90", t4_90, list(
  All         = c(1.79, 1.68, 1.92, 17.9, 15.7, 20.0),
  Transient   = c(2.27, 1.98, 2.61, 17.9, 14.6, 21.3),
  Persistent2 = c(1.97, 1.74, 2.23, 19.8, 16.0, 23.6),
  Worsening   = c(1.36, 1.14, 1.63, 14.5, 6.2, 22.7),
  Severe      = c(1.43, 1.27, 1.60, 15.6, 10.4, 20.7)))
chk_eff("Table 4 hosp", t4_h, list(
  All         = c(1.39, 1.26, 1.54, 5.4, 3.6, 7.2),
  Transient   = c(1.74, 1.39, 2.18, 5.0, 2.7, 7.3),
  Persistent2 = c(1.59, 1.32, 1.91, 7.2, 4.1, 10.2),
  Worsening   = c(1.02, 0.81, 1.30, 1.4, -6.2, 9.0),
  Severe      = c(1.22, 1.03, 1.45, 5.1, 0.3, 9.8)))
RES$table4_d90 <- t4_90; RES$table4_hosp <- t4_h

it90 <- interaction_tests(md3, "d90", "traj")
ith  <- interaction_tests(md3, "hosp_death", "traj")
chk_int("交互 d90", it90, 38.5795, 3.5606)
chk_int("交互 hosp", ith, 15.9265, 3.3122)
RES$interaction_traj <- rbind(cbind(outcome = "d90", it90), cbind(outcome = "hosp", ith))

# 各队列内、相对一过性的 RR / RD(稿子"生存分析"一节与结果第二节)
wc <- md3[md3$cancer == 1, ]; wk <- md3[md3$cancer == 0, ]
tt <- c("trajPersistent2", "trajWorsening", "trajSevere")
within_eff <- function(dat, y, cl = NULL)
  list(rr = tidy_fit(fit_rr(dat, y, c("traj", cov_names), cl), tt),
       rd = tidy_fit(fit_rd(dat, y, c("traj", cov_names), cl), tt))
w1 <- within_eff(wc, "d90")
chk_terms("队列内 肿瘤 d90", w1$rr, w1$rd, list(
  trajPersistent2 = c(1.20, 1.08, 1.34, 6.7, 2.6, 10.8),
  trajWorsening   = c(1.42, 1.23, 1.63, 14.8, 8.2, 21.5),
  trajSevere      = c(1.30, 1.16, 1.46, 11.3, 6.3, 16.3)))
w2 <- within_eff(wk, "d90", "subclass")
chk_terms("队列内 非肿瘤 d90", w2$rr, w2$rd, list(
  trajPersistent2 = c(1.23, 1.08, 1.41, 2.9, 0.7, 5.2),
  trajWorsening   = c(1.97, 1.69, 2.30, 16.7, 12.0, 21.4),
  trajSevere      = c(1.59, 1.39, 1.81, 10.7, 7.5, 14.0)))
w3 <- within_eff(wc, "hosp_death")
chk_terms("队列内 肿瘤 hosp", w3$rr, w3$rd, list(
  trajPersistent2 = c(1.43, 1.17, 1.75, 5.7, 2.4, 8.9),
  trajWorsening   = c(1.94, 1.53, 2.46, 14.6, 8.6, 20.6),
  trajSevere      = c(1.56, 1.26, 1.93, 8.5, 4.2, 12.9)))
w4 <- within_eff(wk, "hosp_death", "subclass")
chk_terms("队列内 非肿瘤 hosp", w4$rr, w4$rd, list(
  trajPersistent2 = c(1.51, 1.24, 1.84, 3.2, 1.4, 5.0),
  trajWorsening   = c(3.10, 2.49, 3.84, 18.4, 13.9, 22.9),
  trajSevere      = c(2.10, 1.73, 2.55, 9.9, 7.2, 12.7)))

# 粗 90 天死亡率(Figure 5)
cr <- crude_risk(md3, "d90", "traj")
exp_cr <- list(Transient = c(14.0, 32.1), Persistent2 = c(19.8, 41.6),
               Worsening = c(39.0, 54.9), Severe = c(35.6, 54.1))
for (s in names(exp_cr)) {
  chk(paste("Figure 5 粗d90% 非肿瘤", s), 100 * cr$risk[cr$stratum == s & cr$exposure == 0], exp_cr[[s]][1], tol = 0.051)
  chk(paste("Figure 5 粗d90% 肿瘤",   s), 100 * cr$risk[cr$stratum == s & cr$exposure == 1], exp_cr[[s]][2], tol = 0.051)
}

# 补充图 S1:加回"早期出 ICU"类 ---------------------------------------------- [V4]
md_alive <- md |>
  select(stay_id, cancer, subclass, all_of(cov_names), hosp_death, d90) |>
  left_join(traj_all |> select(stay_id, pheno6), by = "stay_id") |>
  filter(pheno6 != "1_early_death")
md_alive$traj6 <- factor(ifelse(md_alive$pheno6 %in% four, md_alive$pheno6, "EarlyICUexit"),
                         levels = c(four, "EarlyICUexit"), labels = c(tlab, "EarlyICUexit"))
s1n <- table(md_alive$cancer, md_alive$traj6)
chk("S1 早期出ICU n 对照", s1n["0", "EarlyICUexit"], 2644, tol = 0)
chk("S1 早期出ICU n 肿瘤", s1n["1", "EarlyICUexit"], 1157, tol = 0)
chk("S1 早期出ICU n 合计", sum(s1n[, "EarlyICUexit"]), 3801, tol = 0)
chk("S1 存活至72h n", nrow(md_alive), 11771, tol = 0)
s1e <- cancer_effect(md_alive, "d90", "traj6")
chk_eff("S1 d90", s1e, list(
  All          = c(2.04, 1.92, 2.18, 18.6, 16.9, 20.4),
  EarlyICUexit = c(3.52, 3.01, 4.12, 20.8, 17.9, 23.7)))
chk_int("S1 交互", interaction_tests(md_alive, "d90", "traj6"), 114.7768, 5.1772)
RES$figS1 <- s1e

# 未匹配队列 ---------------------------------------------------------------- [V4]
build_cov <- function(base, outc, cancer) {
  b <- as.data.frame(base); o <- as.data.frame(outc)
  i <- match(b$stay_id, o$stay_id)
  data.frame(stay_id = b$stay_id, cancer = cancer,
    age = b$admission_age, male = as.numeric(b$gender == "M"), sofa = b$sofa_day1,
    sepsis = z0(num01(b$sepsis3)),
    chf   = z0(num01(b$congestive_heart_failure)),
    ckd   = z0(num01(b$renal_disease)),
    copd  = z0(num01(b$chronic_pulmonary_disease)),
    dm    = pmax(z0(num01(b$diabetes_without_cc)), z0(num01(b$diabetes_with_cc))),
    liver = pmax(z0(num01(b$mild_liver_disease)), z0(num01(b$severe_liver_disease))),
    hosp_death = num01(o$hospital_expire_flag[i]),
    d90        = num01(o$death_90d[i]))
}
allp <- rbind(build_cov(baseline_c, outc_c, 1), build_cov(baseline_k, outc_k, 0))
allp <- allp |> left_join(ex_cov, by = "stay_id")                 # v7:新增协变量
stopifnot(!anyNA(allp$adm_type))
allp$pheno6 <- traj_all$pheno6[match(allp$stay_id, traj_all$stay_id)]
# 核对:已匹配的 7970 人,协变量应与 md3 完全一致
mi <- match(md3$stay_id, allp$stay_id)
for (v in c(cov_names, "hosp_death", "d90"))
  chk(paste("未匹配队列 协变量一致", v), as.numeric(isTRUE(all.equal(allp[[v]][mi], md3[[v]]))), 1, tol = 0)
um <- allp[allp$pheno6 %in% four, ]
um <- um[complete.cases(um[, c("age", "male", "sofa")]), ]
um$traj <- factor(um$pheno6, levels = four, labels = tlab)
chk("未匹配 n", nrow(um), 15104, tol = 0)
chk("未匹配 非肿瘤 n", sum(um$cancer == 0), 12448, tol = 0)
chk("未匹配 肿瘤 n", sum(um$cancer == 1), 2656, tol = 0)
# 连续性核对:旧口径 OR 2.92 (2.64–3.23)
f_um <- glm(as.formula(paste("d90 ~ cancer + traj +", covs)), data = um, family = binomial)
ct_um <- coeftest(f_um, vcov. = vcovHC(f_um, type = "HC1"))
chk("未匹配 d90 OR(旧口径)", exp(ct_um["cancer", 1]), 2.94)
chk("未匹配 d90 OR lo", exp(ct_um["cancer", 1] - 1.96 * ct_um["cancer", 2]), 2.63)
chk("未匹配 d90 OR hi", exp(ct_um["cancer", 1] + 1.96 * ct_um["cancer", 2]), 3.28)
um90 <- cancer_effect(um, "d90", "traj", cluster = NULL)
umh  <- cancer_effect(um, "hosp_death", "traj", cluster = NULL)
chk_eff("未匹配 d90", um90, list(
  All         = c(1.84, 1.74, 1.95, 18.3, 16.3, 20.2),
  Transient   = c(2.27, 2.02, 2.55, 17.1, 14.1, 20.2),
  Persistent2 = c(2.00, 1.80, 2.21, 19.7, 16.2, 23.2),
  Worsening   = c(1.52, 1.30, 1.78, 18.3, 11.1, 25.6),
  Severe      = c(1.51, 1.36, 1.67, 17.4, 12.8, 22.1)))
chk_eff("未匹配 hosp", umh, list(All = c(1.45, 1.33, 1.58, 6.0, 4.4, 7.6)))
chk_int("未匹配 交互 d90", interaction_tests(um, "d90", "traj", cluster = NULL), 43.4008, 3.5936)
RES$unmatched_d90 <- um90

# ---- 6. Table 5:按肿瘤类型 --------------------------------------------------- [V4]
sub_c <- baseline_c |>
  transmute(stay_id,
            ctype = case_when(heme == 1 ~ "Hematologic",
                              metastatic_solid_tumor == 1 ~ "Metastatic solid",
                              TRUE ~ "Non-metastatic solid"))
# 对照继承其匹配肿瘤患者的类型(必须从匹配后全部 12,636 人里取,而不是主样本)
sub_map <- md2 |> filter(cancer == 1) |>
  left_join(sub_c, by = "stay_id") |> distinct(subclass, ctype)
md4 <- md3 |> left_join(sub_map, by = "subclass")
stopifnot(!anyNA(md4$ctype))

n5 <- table(md4$cancer, md4$ctype)
exp_n5 <- list(Hematologic = c(979, 514), `Non-metastatic solid` = c(2041, 1041),
               `Metastatic solid` = c(2294, 1101))
for (ct in names(exp_n5)) {
  chk(paste("Table 5 n 对照", ct), n5["0", ct], exp_n5[[ct]][1], tol = 0)
  chk(paste("Table 5 n 肿瘤", ct), n5["1", ct], exp_n5[[ct]][2], tol = 0)
}
t5_90 <- cancer_effect(md4, "d90", "ctype")
t5_h  <- cancer_effect(md4, "hosp_death", "ctype")
chk_eff("Table 5 d90", t5_90, list(
  Hematologic            = c(1.68, 1.42, 1.98, 15.4, 10.3, 20.5),
  `Non-metastatic solid` = c(1.30, 1.15, 1.47, 6.9, 3.5, 10.3),
  `Metastatic solid`     = c(2.39, 2.17, 2.64, 30.2, 26.8, 33.7)))
chk_eff("Table 5 hosp", t5_h, list(
  Hematologic            = c(1.61, 1.26, 2.05, 8.6, 4.2, 13.0),
  `Non-metastatic solid` = c(0.98, 0.82, 1.18, -0.4, -3.1, 2.3),
  `Metastatic solid`     = c(1.74, 1.50, 2.03, 9.6, 6.7, 12.4)))
crude5 <- md4 |> group_by(ctype, cancer) |>
  summarise(hosp = 100 * mean(hosp_death), d90 = 100 * mean(d90), .groups = "drop")
ex5 <- list(Hematologic = c(12.3, 27.0, 20.2, 43.8), `Non-metastatic solid` = c(14.7, 13.5, 23.0, 28.3),
            `Metastatic solid` = c(12.7, 24.2, 21.2, 53.8))
for (ct in names(ex5)) {
  a <- crude5[crude5$ctype == ct & crude5$cancer == 0, ]; b <- crude5[crude5$ctype == ct & crude5$cancer == 1, ]
  chk(paste("Table 5 院内% 非肿瘤", ct), a$hosp, ex5[[ct]][1], tol = 0.051)
  chk(paste("Table 5 院内% 肿瘤",   ct), b$hosp, ex5[[ct]][2], tol = 0.051)
  chk(paste("Table 5 90天% 非肿瘤", ct), a$d90,  ex5[[ct]][3], tol = 0.051)
  chk(paste("Table 5 90天% 肿瘤",   ct), b$d90,  ex5[[ct]][4], tol = 0.051)
}
chk_int("Table 5 异质性 d90", interaction_tests(md4, "d90", "ctype"), 71.1027, 113.0998)
chk_int("Table 5 异质性 hosp", interaction_tests(md4, "hosp_death", "ctype"), 31.0692, 38.0815)
RES$table5_d90 <- t5_90; RES$table5_hosp <- t5_h

# ---- 7. 肾脏结局;Table 6(仅肿瘤患者)与补充表 S1(合并匹配样本)---------------- [V4]
d0 <- bind_rows(daily_c |> select(stay_id, day, max_stage),
                daily_k |> select(stay_id, day, max_stage)) |>
  filter(day == 0) |> transmute(stay_id, day0 = max_stage)

grp <- md3 |> left_join(d0, by = "stay_id") |>
  mutate(g = case_when(
    pheno6 == "3_transient" & day0 < 3  ~ "A_2_transient",
    pheno6 == "3_transient" & day0 >= 3 ~ "B_3_resolved",
    pheno6 == "4_persistent_stage2"     ~ "C_2_persistent",
    pheno6 == "5_worsening"             ~ "D_3_worsening",
    pheno6 == "6_severe_from_start"     ~ "E_3_severe"))
stopifnot(!anyNA(grp$g))

kid <- bind_rows(kid_fun(outc_c, baseline_c), kid_fun(outc_k, baseline_k))
grp2 <- grp |> left_join(kid, by = "stay_id")
grp2$g5     <- relevel(factor(grp2$g), ref = "A_2_transient")   # 参照:峰值 2 期、一过性
grp2$peak3  <- as.numeric(grp2$g %in% c("B_3_resolved", "D_3_worsening", "E_3_severe"))
grp2$unres  <- as.numeric(grp2$g %in% c("C_2_persistent", "D_3_worsening", "E_3_severe"))
G <- c("A_2_transient", "B_3_resolved", "C_2_persistent", "D_3_worsening", "E_3_severe")

summ_g <- function(d) {
  s <- d |> group_by(g) |>
    summarise(n = n(), d90 = 100 * mean(d90),
              rrt_late = 100 * mean(y_rrt_late, na.rm = TRUE),
              nonrec = 100 * mean(y_nonrec, na.rm = TRUE),
              make = 100 * mean(y_make), .groups = "drop") |> as.data.frame()
  rownames(s) <- s$g; s
}
# Table 6:仅肿瘤患者(n = 2,656)的各组百分比
ca <- grp2[grp2$cancer == 1, ]
s6c <- summ_g(ca)
for (i in seq_along(G)) {
  chk(paste("Table 6 n", G[i]),        s6c[G[i], "n"],        c(891, 98, 866, 235, 566)[i], tol = 0)
  chk(paste("Table 6 d90%", G[i]),     s6c[G[i], "d90"],      c(31.0, 41.8, 41.6, 54.9, 54.1)[i], tol = 0.051)
  chk(paste("Table 6 72h后RRT%", G[i]), s6c[G[i], "rrt_late"], c(0.9, 2.1, 1.7, 4.6, 5.3)[i], tol = 0.051)
  chk(paste("Table 6 非恢复%", G[i]),  s6c[G[i], "nonrec"],   c(6.9, 9.2, 12.0, 23.6, 31.0)[i], tol = 0.051)
  chk(paste("Table 6 MAKE%", G[i]),    s6c[G[i], "make"],     c(15.0, 18.4, 24.8, 44.3, 47.5)[i], tol = 0.051)
}
# Table 6:相对 A_2_transient 的 RR / RD(观测独立,不聚类)
e6 <- list(
  d90 = list(B_3_resolved   = c(1.21, 0.96, 1.52, 6.6, -2.6, 15.7),
             C_2_persistent = c(1.23, 1.09, 1.38, 7.3, 3.2, 11.5),
             D_3_worsening  = c(1.45, 1.26, 1.68, 15.5, 8.8, 22.2),
             E_3_severe     = c(1.33, 1.18, 1.51, 12.0, 6.9, 17.1)),
  y_make = list(B_3_resolved   = c(1.17, 0.77, 1.77, 1.9, -5.7, 9.4),
                C_2_persistent = c(1.51, 1.25, 1.82, 8.0, 4.4, 11.5),
                D_3_worsening  = c(2.19, 1.78, 2.68, 21.6, 15.5, 27.8),
                E_3_severe     = c(2.10, 1.75, 2.52, 20.6, 15.8, 25.5)),
  y_nonrec = list(B_3_resolved   = c(1.47, 0.70, 3.07, 3.5, -3.1, 10.2),
                  C_2_persistent = c(1.72, 1.22, 2.42, 5.1, 1.9, 8.3),
                  D_3_worsening  = c(3.05, 2.04, 4.55, 15.4, 8.5, 22.3),
                  E_3_severe     = c(3.35, 2.41, 4.66, 19.3, 14.1, 24.5)),
  y_rrt_late = list(B_3_resolved   = c(2.40, 0.50, 11.57, 1.0, -2.0, 4.0),
                    C_2_persistent = c(1.99, 0.85, 4.70, 0.9, -0.1, 2.0),
                    D_3_worsening  = c(3.56, 1.39, 9.12, 3.1, 0.2, 5.9),
                    E_3_severe     = c(3.95, 1.78, 8.78, 3.1, 1.1, 5.0)))
RES$table6_effects <- list()
for (y in names(e6)) {
  ge <- group_effects(ca, y, "g5")
  chk_grp(paste("Table 6", y), ge, e6[[y]])
  RES$table6_effects[[y]] <- ge
}

# 补充表 S1:合并匹配样本(n = 7,970),调整肿瘤状态,按匹配集聚类
pooled_effects <- function(data, y) {
  tt5 <- paste0("g5", levels(data$g5)[-1])
  rbind(cbind(scale = "RR", tidy_fit(fit_rr(data, y, c("g5", "cancer", cov_names), "subclass"), tt5)),
        cbind(scale = "RD", tidy_fit(fit_rd(data, y, c("g5", "cancer", cov_names), "subclass"), tt5)))
}
s6p <- summ_g(grp2)
for (i in seq_along(G)) {
  chk(paste("S1 n", G[i]),        s6p[G[i], "n"],        c(2907, 276, 2626, 620, 1541)[i], tol = 0)
  chk(paste("S1 d90%", G[i]),     s6p[G[i], "d90"],      c(18.9, 27.9, 27.0, 45.0, 42.4)[i], tol = 0.051)
  chk(paste("S1 72h后RRT%", G[i]), s6p[G[i], "rrt_late"], c(0.6, 2.3, 1.6, 8.4, 6.7)[i], tol = 0.051)
  chk(paste("S1 非恢复%", G[i]),  s6p[G[i], "nonrec"],   c(6.8, 8.8, 11.1, 23.8, 29.1)[i], tol = 0.051)
  chk(paste("S1 MAKE%", G[i]),    s6p[G[i], "make"],     c(10.0, 17.0, 19.0, 43.7, 44.1)[i], tol = 0.051)
}
ep <- list(
  d90 = list(B_3_resolved   = c(1.19, 0.99, 1.43, 2.7, -2.2, 7.7),
             C_2_persistent = c(1.26, 1.15, 1.38, 4.6, 2.5, 6.7),
             D_3_worsening  = c(1.72, 1.54, 1.93, 16.1, 12.2, 20.0),
             E_3_severe     = c(1.49, 1.35, 1.63, 11.3, 8.5, 14.0)),
  y_make = list(B_3_resolved   = c(1.47, 1.12, 1.91, 3.7, -0.5, 7.9),
                C_2_persistent = c(1.72, 1.51, 1.96, 6.9, 5.1, 8.8),
                D_3_worsening  = c(3.13, 2.73, 3.60, 25.9, 22.1, 29.8),
                E_3_severe     = c(2.78, 2.45, 3.16, 22.9, 20.3, 25.6)),
  y_nonrec = list(B_3_resolved   = c(1.43, 0.93, 2.21, 3.0, -0.7, 6.7),
                  C_2_persistent = c(1.62, 1.34, 1.95, 4.2, 2.5, 5.9),
                  D_3_worsening  = c(3.03, 2.43, 3.78, 15.2, 11.0, 19.3),
                  E_3_severe     = c(3.24, 2.71, 3.87, 18.1, 15.2, 21.0)),
  y_rrt_late = list(B_3_resolved   = c(3.53, 1.46, 8.55, 1.5, -0.3, 3.2),
                    C_2_persistent = c(2.58, 1.47, 4.53, 0.9, 0.3, 1.4),
                    D_3_worsening  = c(9.38, 5.37, 16.41, 6.7, 4.5, 8.9),
                    E_3_severe     = c(6.62, 3.91, 11.22, 4.4, 3.1, 5.6)))
RES$tableS1_effects <- list()
for (y in names(ep)) {
  pe_ <- pooled_effects(grp2, y)
  chk_grp(paste("S1", y), pe_, ep[[y]])
  RES$tableS1_effects[[y]] <- pe_
}
chk_int("S1 cancer×组 d90",  interaction_tests(grp2, "d90", "g5"), 38.3138, 4.8175)
chk_int("S1 cancer×组 make", interaction_tests(grp2, "y_make", "g5"), 27.3764, 4.7026)

# ---- 7b. 峰值分期之外,轨迹是否带来额外信息(仅肿瘤患者)------------------------ [V4]
m0 <- glm(as.formula(paste("d90 ~", covs)), data = ca, family = binomial)
m1 <- glm(as.formula(paste("d90 ~ peak3 +", covs)), data = ca, family = binomial)
m2 <- glm(as.formula(paste("d90 ~ factor(g) +", covs)), data = ca, family = binomial)
lrt <- anova(m1, m2, test = "Chisq")
chk("LRT 轨迹加入峰值分期 chisq", lrt$Deviance[2], 14.816, tol = 0.01)
chk("LRT 轨迹加入峰值分期 df", lrt$Df[2], 3, tol = 0)
chk("LRT 轨迹加入峰值分期 P", lrt[["Pr(>Chi)"]][2], 0.0020, tol = 0.00051)

# 峰值分期内:未缓解 vs 缓解,及其与峰值分期的交互
ur_rr <- tidy_fit(fit_rr(ca, "d90", c("peak3*unres", cov_names)), c("unres", "peak3:unres"))
ur_rd <- tidy_fit(fit_rd(ca, "d90", c("peak3*unres", cov_names)), c("unres", "peak3:unres"))
chk_terms("峰值2期内 持续 vs 一过性", ur_rr, ur_rd, list(
  unres         = c(1.23, 1.09, 1.38, 7.4, 3.2, 11.5),
  `peak3:unres` = c(0.92, 0.72, 1.18, -0.9, -11.0, 9.3)))
chk("峰值×未缓解 交互P(RR)", ur_rr$p[2], 0.52, tol = 0.0051)
chk("峰值×未缓解 交互P(RD)", ur_rd$p[2], 0.87, tol = 0.0051)
s3 <- ca[ca$peak3 == 1, ]
chk_vec("峰值3期内 未缓解 vs 缓解", c(1.18, 0.95, 1.47, 6.9, -2.7, 16.5),
        tidy_fit(fit_rr(s3, "d90", c("unres", cov_names)), "unres"),
        tidy_fit(fit_rd(s3, "d90", c("unres", cov_names)), "unres"))

# 五折交叉验证 AUC:折在合并匹配样本(7,970 例)上按匹配集随机分配(set.seed(1)),
# 再取肿瘤患者;与 v3 稿子及 section15.R 的做法一致。AUC 及差值的 bootstrap 区间
# 依赖这一折分配,改变分折方式会使区间下限在 0 附近波动。
set.seed(1)
ss <- unique(grp2$subclass)
fold_of <- setNames(sample(rep(1:5, length.out = length(ss))), ss)
ca$fold <- fold_of[as.character(ca$subclass)]
forms <- list(m0 = paste("d90 ~", covs),
              m1 = paste("d90 ~ peak3 +", covs),
              m2 = paste("d90 ~ factor(g) +", covs))
cvpred <- function(f) {
  p <- numeric(nrow(ca))
  for (k in 1:5) {
    fit <- glm(as.formula(f), data = ca[ca$fold != k, ], family = binomial)
    p[ca$fold == k] <- predict(fit, newdata = ca[ca$fold == k, ], type = "response")
  }
  p
}
preds <- lapply(forms, cvpred)
aucs  <- sapply(preds, function(p) as.numeric(auc(roc(ca$d90, p, quiet = TRUE))))
nt_auc <- "折在合并匹配样本上分配(set.seed(1)),bootstrap set.seed(2);与 section15.R(v3)相同"
chk("AUC 仅协变量",       aucs["m0"], 0.735, tol = 0.0006, note = nt_auc)
chk("AUC 加峰值分期",     aucs["m1"], 0.740, tol = 0.0006, note = nt_auc)
chk("AUC 加峰值分期+轨迹", aucs["m2"], 0.743, tol = 0.0006, note = nt_auc)
set.seed(2)
nb <- nrow(ca)
dauc <- replicate(2000, {
  i <- sample(nb, nb, replace = TRUE)
  as.numeric(auc(roc(ca$d90[i], preds$m2[i], quiet = TRUE))) -
    as.numeric(auc(roc(ca$d90[i], preds$m1[i], quiet = TRUE)))
})
ci_d <- quantile(dauc, c(0.025, 0.975))
chk("AUC 差值", aucs["m2"] - aucs["m1"], 0.0030, tol = 0.0006, note = nt_auc)
chk("AUC 差值 bootstrap lo", ci_d[1], -0.0009, tol = 0.0006, note = nt_auc)
chk("AUC 差值 bootstrap hi", ci_d[2], 0.0069, tol = 0.0006, note = nt_auc)
RES$auc <- data.frame(model = names(aucs), auc = round(aucs, 4))

# ---- 8. Table 7:肿瘤 vs 非肿瘤的肾脏结局 ------------------------------------ [V4]
kd <- md3 |> left_join(kid, by = "stay_id")
ex7 <- list(
  any_rrt  = list(y = "y_rrt_any",  e = c(0.94, 0.78, 1.13, -0.7, -1.7, 0.4), n = c(2656, 149, 5314, 331)),
  rrt_late = list(y = "y_rrt_late", e = c(0.81, 0.59, 1.11, -0.6, -1.4, 0.1), n = c(2569,  62, 5124, 141)),
  nonrec   = list(y = "y_nonrec",   e = c(1.20, 1.05, 1.38, 2.4, 0.6, 4.3), n = c(1894, 272, 4162, 522)),
  make     = list(y = "y_make",     e = c(1.29, 1.20, 1.40, 6.0, 4.1, 7.9), n = c(2656, 740, 5314, 1049)))
t7 <- NULL
for (k in names(ex7)) {
  y <- ex7[[k]]$y
  rr <- tidy_fit(fit_rr(kd, y, c("cancer", "traj", cov_names), "subclass"), "cancer")
  rd <- tidy_fit(fit_rd(kd, y, c("cancer", "traj", cov_names), "subclass"), "cancer")
  chk_vec(paste("Table 7", k), ex7[[k]]$e, rr, rd)
  d <- kd[!is.na(kd[[y]]), ]
  chk(paste("Table 7", k, "n 肿瘤"),     sum(d$cancer == 1), ex7[[k]]$n[1], tol = 0)
  chk(paste("Table 7", k, "事件 肿瘤"),   sum(d[[y]][d$cancer == 1]), ex7[[k]]$n[2], tol = 0)
  chk(paste("Table 7", k, "n 非肿瘤"),   sum(d$cancer == 0), ex7[[k]]$n[3], tol = 0)
  chk(paste("Table 7", k, "事件 非肿瘤"), sum(d[[y]][d$cancer == 0]), ex7[[k]]$n[4], tol = 0)
  t7 <- rbind(t7, data.frame(outcome = k, RR = rr$estimate, RR_lo = rr$lo, RR_hi = rr$hi,
                             RD_pts = 100 * rd$estimate, RD_lo = 100 * rd$lo, RD_hi = 100 * rd$hi))
}
RES$table7 <- t7

# ---- 10. Table 1 要点(肿瘤队列,各轨迹类的事件数)----------------------------- [R]
t1 <- outc_c |> left_join(traj_c |> select(stay_id, pheno6), by = "stay_id") |>
  group_by(pheno6) |>
  summarise(hosp = sum(to01(hospital_expire_flag)), d28 = sum(to01(death_28d)),
            d90 = sum(to01(death_90d)), rrt = sum(to01(any_rrt)), .groups = "drop") |>
  as.data.frame()
rownames(t1) <- t1$pheno6
cl <- rownames(t1)
ex1 <- list(hosp = c(388, 76, 124, 173, 78, 171), d28 = c(399, 208, 200, 244, 101, 236),
            d90  = c(399, 362, 317, 360, 129, 306), rrt = c(18, 2, 15, 20, 29, 85))
for (v in names(ex1)) for (i in 1:6)
  chk(paste("Table 1", v, cl[i]), t1[cl[i], v], ex1[[v]][i], tol = 0)
RES$table1_counts <- t1

# ---- 11. RRT 是否导致轻度轨迹类的误分类 ---------------------------------------- [R]
rrt_all <- bind_rows(outc_c |> transmute(stay_id, rrt_h = as.numeric(rrt_h)),
                     outc_k |> transmute(stay_id, rrt_h = as.numeric(rrt_h)))
rr <- md2 |> left_join(rrt_all, by = "stay_id") |>
  mutate(light = pheno6 %in% c("3_transient", "4_persistent_stage2"),
         rrt72 = !is.na(rrt_h) & rrt_h <= 72) |>
  group_by(cancer) |>
  summarise(n_light = sum(light), n_light_rrt72 = sum(light & rrt72),
            pct = round(100 * n_light_rrt72 / n_light, 2), n_rrt72_all = sum(rrt72),
            .groups = "drop") |> as.data.frame()
chk("RRT核查 肿瘤 轻类n", rr$n_light[rr$cancer == 1], 1855, tol = 0)
chk("RRT核查 肿瘤 轻类中前72hRRT", rr$n_light_rrt72[rr$cancer == 1], 10, tol = 0)
chk("RRT核查 非肿瘤 轻类n", rr$n_light[rr$cancer == 0], 3954, tol = 0)
chk("RRT核查 非肿瘤 轻类中前72hRRT", rr$n_light_rrt72[rr$cancer == 0], 14, tol = 0)
chk("前72h开始RRT总数 肿瘤", rr$n_rrt72_all[rr$cancer == 1], 107, tol = 0)
RES$rrt_misclass <- rr

# t0 是否未受 RRT 影响 ------------------------------------------------------------ [R]
t0chk <- kdigo_temp |> filter(hours_from_t0 == 0) |>
  mutate(st = pmax(aki_stage_creat, aki_stage_uo, na.rm = TRUE)) |>
  summarise(n = n(), n_ge2 = sum(st >= 2, na.rm = TRUE), n_crrt = sum(!is.na(aki_stage_crrt)))
chk("t0核查 肿瘤组 肌酐/尿量≥2期 = n", t0chk$n_ge2, t0chk$n, tol = 0)
chk("t0核查 肿瘤组 t0时有CRRT记录", t0chk$n_crrt, 0, tol = 0)
chk("t0核查 非肿瘤 第0天均≥2期",
    sum(daily_k$max_stage[daily_k$day == 0] >= 2, na.rm = TRUE), 19335, tol = 0)

# rrt_h 的取整方式:只用 outc 自身的 rrt_time 与 t0,不依赖 procedureevents ------- [N]
rt <- bind_rows(
  outc_c |> transmute(t0 = to_time(t0), rrt_time = to_time(rrt_time), rrt_h = as.numeric(rrt_h)),
  outc_k |> transmute(t0 = to_time(t0), rrt_time = to_time(rrt_time), rrt_h = as.numeric(rrt_h))) |>
  filter(!is.na(rrt_time), !is.na(rrt_h))
h_raw <- as.numeric(difftime(rt$rrt_time, rt$t0, units = "hours"))
# 稿子写明"按钟点之差计算":先把两个时间截到整点再相减,应与 rrt_h 完全一致
h_clock <- as.numeric(difftime(trunc(rt$rrt_time, "hours"), trunc(rt$t0, "hours"), units = "hours"))
chk("rrt_h = 钟点之差(一致比例)", mean(h_clock == rt$rrt_h), 1, tol = 0)
cat("\n---- rrt_h 取整方式(rrt_time 与 t0 之差,共", nrow(rt), "例)----\n",
    "  round  :", round(mean(round(h_raw) == rt$rrt_h), 4), "\n",
    "  floor  :", round(mean(floor(h_raw) == rt$rrt_h), 4), "\n",
    "  ceiling:", round(mean(ceiling(h_raw) == rt$rrt_h), 4), "\n", sep = "")

# ---- 12. 出院去向(安宁疗护)敏感性分析 ---------------------------------------- [V4]
find_mimic <- function(explicit, stem) {
  if (!is.na(explicit) && file.exists(explicit)) return(explicit)
  if (is.na(MIMIC_DIR) || !dir.exists(MIMIC_DIR)) return(NA_character_)
  hit <- list.files(MIMIC_DIR, pattern = paste0("^", stem, "\\.csv(\\.gz)?$"),
                    recursive = TRUE, full.names = TRUE)
  if (length(hit) > 0) hit[1] else NA_character_
}
adm_path <- find_mimic(ADM_PATH, "admissions")
pe_path  <- find_mimic(PE_PATH,  "procedureevents")
di_path  <- find_mimic(DI_PATH,  "d_items")
cat("\n---- MIMIC 文件查找结果 ----\n",
    "admissions:      ", ifelse(is.na(adm_path), "未找到", adm_path), "\n",
    "procedureevents: ", ifelse(is.na(pe_path),  "未找到", pe_path),  "\n",
    "d_items:         ", ifelse(is.na(di_path),  "未找到", di_path),  "\n", sep = "")
if (is.na(adm_path) || is.na(pe_path))
  warning("第 12、13 节需要 admissions 和 procedureevents;请设置 MIMIC_DIR 或 ADM_PATH / PE_PATH")
if (!is.na(adm_path)) {
  adm <- read_gz_cols(adm_path, c("hadm_id", "discharge_location"))
  adm$hadm_id <- as.numeric(adm$hadm_id)
  cat("\n---- discharge_location 取值(请确认安宁疗护对应哪个标签)----\n")
  print(sort(table(adm$discharge_location, useNA = "ifany"), decreasing = TRUE))
  adm$hospice <- grepl("HOSPICE", toupper(adm$discharge_location))
  hmap <- bind_rows(outc_c |> transmute(stay_id, hadm_id = as.numeric(hadm_id)),
                    outc_k |> transmute(stay_id, hadm_id = as.numeric(hadm_id)))
  kd2 <- kd |> left_join(hmap, by = "stay_id") |>
    left_join(adm |> select(hadm_id, hospice) |> distinct(hadm_id, .keep_all = TRUE), by = "hadm_id") |>
    mutate(hospice = hospice %in% TRUE)
  surv <- kd2 |> filter(hosp_death == 0) |> group_by(cancer) |>
    summarise(n = n(), n_hospice = sum(hospice), pct = 100 * mean(hospice), .groups = "drop")
  chk("安宁疗护% 存活者 肿瘤", surv$pct[surv$cancer == 1], 11.3, tol = 0.06)
  chk("安宁疗护% 存活者 非肿瘤", surv$pct[surv$cancer == 0], 2.9, tol = 0.06)
  chk("排除安宁疗护的人数", sum(kd2$hospice), 372, tol = 0)
  kdn <- kd2 |> filter(!hospice)
  chk("排除后 n", nrow(kdn), 7598, tol = 0)
  chk("排除后 肿瘤 n", sum(kdn$cancer == 1), 2416, tol = 0)
  chk("排除后 非肿瘤 n", sum(kdn$cancer == 0), 5182, tol = 0)
  hs <- list(
    d90    = list(y = "d90",      e = c(1.75, 1.62, 1.88, 15.6, 13.3, 17.8), n = c(2416, 900, 5182, 1039)),
    make   = list(y = "y_make",   e = c(1.32, 1.22, 1.43, 6.8, 4.8, 8.8), n = c(2416, 692, 5182, 1030)),
    nonrec = list(y = "y_nonrec", e = c(1.15, 0.99, 1.32, 1.8, -0.1, 3.6), n = c(1722, 227, 4062, 497)))
  s7 <- NULL
  for (k in names(hs)) {
    y <- hs[[k]]$y
    rr <- tidy_fit(fit_rr(kdn, y, c("cancer", "traj", cov_names), "subclass"), "cancer")
    rd <- tidy_fit(fit_rd(kdn, y, c("cancer", "traj", cov_names), "subclass"), "cancer")
    chk_vec(paste("排除安宁疗护", k), hs[[k]]$e, rr, rd)
    d <- kdn[!is.na(kdn[[y]]), ]
    chk(paste("排除安宁疗护", k, "n 肿瘤"),     sum(d$cancer == 1), hs[[k]]$n[1], tol = 0)
    chk(paste("排除安宁疗护", k, "事件 肿瘤"),   sum(d[[y]][d$cancer == 1]), hs[[k]]$n[2], tol = 0)
    chk(paste("排除安宁疗护", k, "n 非肿瘤"),   sum(d$cancer == 0), hs[[k]]$n[3], tol = 0)
    chk(paste("排除安宁疗护", k, "事件 非肿瘤"), sum(d[[y]][d$cancer == 0]), hs[[k]]$n[4], tol = 0)
    s7 <- rbind(s7, data.frame(outcome = k, RR = rr$estimate, RR_lo = rr$lo, RR_hi = rr$hi,
                               RD_pts = 100 * rd$estimate, RD_lo = 100 * rd$lo, RD_hi = 100 * rd$hi))
  }
  RES$hospice_survivors <- as.data.frame(surv); RES$table7_no_hospice <- s7
} else message("【跳过第 12 节】未找到 admissions.csv(.gz),安宁疗护敏感性分析未运行")

# ---- 12b. worsening vs severe-from-start(补充材料)---------------------------- [V4]
# 需要 SOFA_PATH:stay_id, sofa, renal(mimiciv_derived.first_day_sofa 在 BigQuery 的导出)
if (!is.na(SOFA_PATH) && file.exists(SOFA_PATH)) {
  sr <- read.csv(SOFA_PATH)
  ws <- md3[md3$traj %in% c("Worsening", "Severe"), ]
  ws$worse <- as.numeric(ws$traj == "Worsening")
  ws$sofa_renal <- sr$renal[match(ws$stay_id, sr$stay_id)]
  chk("补充 SOFA 与导出文件一致", mean(ws$sofa == sr$sofa[match(ws$stay_id, sr$stay_id)], na.rm = TRUE), 1, tol = 0)
  chk("补充 worsening+severe n", nrow(ws), 2161, tol = 0)
  ws$sofa_nr <- ws$sofa - ws$sofa_renal
  covs_nr <- sub("^sofa$", "sofa_nr", cov_names)
  chk_vec("补充 总SOFA", c(1.16, 1.05, 1.27, 5.7, 1.4, 9.9),
          tidy_fit(fit_rr(ws, "d90", c("worse", "cancer", cov_names), "subclass"), "worse"),
          tidy_fit(fit_rd(ws, "d90", c("worse", "cancer", cov_names), "subclass"), "worse"))
  chk_vec("补充 去肾脏SOFA", c(1.12, 1.02, 1.24, 4.6, 0.3, 8.8),
          tidy_fit(fit_rr(ws, "d90", c("worse", "cancer", covs_nr), "subclass"), "worse"),
          tidy_fit(fit_rd(ws, "d90", c("worse", "cancer", covs_nr), "subclass"), "worse"))
} else message("【跳过第 12b 节】未设置 SOFA_PATH(stay_id, sofa, renal 的 csv)")

# ---- 13. RRT 数据来源核查(procedureevents)------------------------------------ [N]
# 目的:找出哪一组 itemid 能复现 outc 里的 any_rrt(肿瘤组 169 例),
#       以及 rrt_h 是按"四舍五入 / 向下 / 向上"取整。结果用来写方法部分的 RRT 来源。
if (!is.na(pe_path)) {
  outc_all <- bind_rows(
    outc_c |> transmute(stay_id, t0 = to_time(t0), rrt_time = to_time(rrt_time),
                        any_rrt = to01(any_rrt), rrt_h = as.numeric(rrt_h), cohort = "cancer"),
    outc_k |> transmute(stay_id, t0 = to_time(t0), rrt_time = to_time(rrt_time),
                        any_rrt = to01(any_rrt), rrt_h = as.numeric(rrt_h), cohort = "non-cancer"))
  cand <- c(225441, 225802, 225803, 225805, 224270, 225809, 225955)   # 凭记忆的备选,以 d_items 为准
  cath <- 224270; crrt_ids <- c(225802, 225803, 225809, 225955)
  if (!is.na(di_path)) {
    di <- read_gz_cols(di_path, c("itemid", "label", "linksto"))
    di_c <- di[di$linksto == "procedureevents" &
               grepl("dialy|crrt|cvvh|hemofil|scuf|ultrafil", di$label, ignore.case = TRUE), ]
    cat("\n---- d_items 中与透析相关的 procedureevents 条目 ----\n"); print(di_c)
    cand     <- as.numeric(di_c$itemid)
    cath     <- as.numeric(di_c$itemid[grepl("catheter", di_c$label, ignore.case = TRUE)])
    crrt_ids <- as.numeric(di_c$itemid[grepl("crrt|cvvh|scuf", di_c$label, ignore.case = TRUE)])
  }
  pe <- read_gz_cols(pe_path, c("stay_id", "starttime", "itemid"))
  pe$stay_id <- as.numeric(pe$stay_id); pe$itemid <- as.numeric(pe$itemid)
  pe <- pe[pe$itemid %in% cand & pe$stay_id %in% outc_all$stay_id, ]
  pe$starttime <- to_time(pe$starttime)

  first_start <- function(items) {
    pe[pe$itemid %in% items, ] |> group_by(stay_id) |>
      summarise(first_start = suppressWarnings(min(starttime)), .groups = "drop")
  }
  eval_set <- function(items, label) {
    d <- outc_all |> left_join(first_start(items), by = "stay_id")
    has <- !is.na(d$first_start)
    dt  <- abs(as.numeric(difftime(d$first_start, d$rrt_time, units = "mins")))
    data.frame(set = label, n_set = sum(has),
               both = sum(has & d$any_rrt == 1), set_only = sum(has & d$any_rrt == 0),
               any_rrt_only = sum(!has & d$any_rrt == 1),
               same_time_1h = sum(has & d$any_rrt == 1 & !is.na(dt) & dt <= 60))
  }
  sets <- list(`all candidate items` = cand,
               `without dialysis catheter` = setdiff(cand, cath),
               `CRRT items only` = crrt_ids)
  for (it in cand) sets[[paste("item", it)]] <- it
  aud <- do.call(rbind, lapply(names(sets), function(s) eval_set(sets[[s]], s)))
  cat("\n---- any_rrt 与 procedureevents 条目的对应关系(全部 23,548 人;any_rrt = TRUE 共",
      sum(outc_all$any_rrt), "例)----\n")
  print(aud, row.names = FALSE)
  best <- aud[order(aud$any_rrt_only + aud$set_only), ][1, ]
  cat("\n最接近的一组:", best$set, "\n")

  d <- outc_all |> left_join(first_start(sets[[best$set]]), by = "stay_id") |>
    filter(any_rrt == 1, !is.na(first_start), !is.na(rrt_h))
  h <- as.numeric(difftime(d$first_start, d$t0, units = "hours"))
  cat("rrt_h 取整方式(与最接近那组 itemid 的开始时间比较):\n")
  cat("  round  :", mean(round(h) == d$rrt_h), "\n")
  cat("  floor  :", mean(floor(h) == d$rrt_h), "\n")
  cat("  ceiling:", mean(ceiling(h) == d$rrt_h), "\n")
  RES$rrt_source_audit <- aud
} else message("【跳过第 13 节】未找到 procedureevents.csv(.gz),RRT 来源核查未运行")

# ---- 15. v7 补充:无条件估计、逐步缩小样本、早期出 ICU 拆分、校正集 A/B/C、E 值 ------------- [V7]
# A = 仅 9 个匹配协变量;B = A + 入院类型、首个 ICU 类型、入院科室、年份段(主分析);C = B + t0 前 code status(探索性)
adj_sets <- list(A = cov_match, B = cov_names, C = c(cov_names, "code"))
fit_est <- function(data, y, adj, cl = "subclass") {
  rr <- tidy_fit(fit_rr(data, y, c("cancer", adj), cl), "cancer")
  rd <- tidy_fit(fit_rd(data, y, c("cancer", adj), cl), "cancer")
  data.frame(RR = rr$estimate, RR_lo = rr$lo, RR_hi = rr$hi,
             RD_pts = 100 * rd$estimate, RD_lo = 100 * rd$lo, RD_hi = 100 * rd$hi)
}
evalue <- function(rr) { rr <- ifelse(rr < 1, 1 / rr, rr); rr + sqrt(rr * (rr - 1)) }   # VanderWeele & Ding 2017

# 15a. 全部匹配患者(12,636 人;72 h 内死亡计为 90 天死亡;不校正轨迹)
stopifnot(nrow(md) == 12636)
v7_uncond <- do.call(rbind, lapply(c("d90", "hosp_death"), function(y)
  do.call(rbind, lapply(names(adj_sets), function(s) cbind(outcome = y, set = s, fit_est(md, y, adj_sets[[s]]))))))
v7_uncond$E_point <- round(evalue(v7_uncond$RR), 2)
v7_uncond$E_limit <- ifelse(v7_uncond$RR_lo > 1, round(evalue(v7_uncond$RR_lo), 2), 1)
RES$v7_uncond <- v7_uncond
RES$v7_crude_all <- aggregate(cbind(d90, hosp_death) ~ cancer, data = md, FUN = mean)

# 15b. 逐步缩小样本(校正集 B):全部 → 排除早期死亡 → 排除早期出 ICU → 校正轨迹
v7_df <- md |> left_join(traj_all |> select(stay_id, pheno6), by = "stay_id")
v7_steps <- list(
  "A all matched"                = v7_df,
  "B excl early death"           = filter(v7_df, pheno6 != "1_early_death"),
  "C excl early ICU exit only"   = filter(v7_df, pheno6 != "2_early_icu_exit"),
  "D excl both (trajectory sample)" = filter(v7_df, pheno6 %in% four))
v7_st <- do.call(rbind, lapply(names(v7_steps), function(nm)
  cbind(sample = nm, n = nrow(v7_steps[[nm]]), fit_est(v7_steps[[nm]], "d90", cov_names))))
v7_st <- rbind(v7_st, cbind(sample = "E trajectory sample + trajectory adj", n = nrow(v7_steps[[4]]),
                            fit_est(v7_steps[[4]], "d90", c(cov_names, "pheno6"))))
RES$v7_steps <- v7_st

# 15c. 早期死亡但出院时 hosp_death = 0 的人数(按癌症状态)
v7_ed <- v7_df[v7_df$pheno6 == "1_early_death", ]
RES$v7_early_death_alive_at_discharge <- as.data.frame(table(cancer = v7_ed$cancer, hosp_death = v7_ed$hosp_death))

# 15d. 早期出 ICU 类按 t0 到转出 ICU 的时间拆分(<24 h;24–72 h;仍在 ICU 但无数据)
oc7 <- bind_rows(
  transmute(outc_c, stay_id, icu_h = as.numeric(difftime(to_time(icu_outtime), to_time(t0), units = "hours"))),
  transmute(outc_k, stay_id, icu_h = as.numeric(difftime(to_time(icu_outtime), to_time(t0), units = "hours"))))
v7_ee <- v7_df |> filter(pheno6 == "2_early_icu_exit") |> left_join(oc7, by = "stay_id") |>
  mutate(sub = case_when(is.na(icu_h) ~ "2x_unknown", icu_h < 24 ~ "2a_left_ICU_<24h",
                         icu_h < 72 ~ "2b_left_ICU_24-72h", TRUE ~ "2c_in_ICU_>=72h_no_data"))
RES$v7_ee_n <- as.data.frame(table(sub = v7_ee$sub, cancer = v7_ee$cancer))
RES$v7_ee_crude <- v7_ee |> group_by(sub, cancer) |>
  summarise(n = n(), hosp = 100 * mean(hosp_death), d90 = 100 * mean(d90), .groups = "drop")
v7_ee_fit <- do.call(rbind, lapply(sort(unique(v7_ee$sub)), function(s) {
  d <- v7_ee[v7_ee$sub == s, ]
  tryCatch(suppressWarnings(cbind(sub = s, n = nrow(d), fit_est(d, "d90", cov_names))),   # 2c 仅 6 人,方差不可估计
           error = function(e) data.frame(sub = s, n = nrow(d), RR = NA, RR_lo = NA, RR_hi = NA,
                                          RD_pts = NA, RD_lo = NA, RD_hi = NA))
}))
v7_ee_fit <- rbind(v7_ee_fit, cbind(sub = "all early ICU exit", n = nrow(v7_ee), fit_est(v7_ee, "d90", cov_names)))
RES$v7_ee_fit <- v7_ee_fit
v7_base <- v7_df |> filter(pheno6 %in% four) |> mutate(grp = pheno6)
sens7 <- function(keep, label) {
  d <- bind_rows(v7_base, mutate(filter(v7_ee, sub %in% keep), grp = "2_early_icu_exit"))
  cbind(added = label, n = nrow(d), fit_est(d, "d90", c(cov_names, "grp")))
}
RES$v7_ee_sens <- rbind(
  sens7(character(0), "none (trajectory sample)"),
  sens7(c("2a_left_ICU_<24h", "2b_left_ICU_24-72h"), "2a + 2b (left ICU)"),
  sens7("2c_in_ICU_>=72h_no_data", "2c only (in ICU, no data)"),
  sens7(unique(v7_ee$sub), "all early ICU exit"))

# 15e. 各轨迹内 90 天 RD / RR,三种校正集(轨迹样本)
v7_tr <- do.call(rbind, lapply(names(adj_sets), function(s)
  cbind(set = s, cancer_effect(v7_base, "d90", "pheno6", adj = adj_sets[[s]]))))
RES$v7_traj_sets <- v7_tr

# 15f. 新协变量的平衡(匹配后全部 12,636 人)与 code status 记录情况
smd_bin <- function(x, g) { x <- as.numeric(x); s <- sqrt((var(x[g == 1]) + var(x[g == 0])) / 2); (mean(x[g == 1]) - mean(x[g == 0])) / s }
RES$v7_balance <- do.call(rbind, lapply(c("adm_type", "unit", "svc", "period", "code"), function(v) {
  lv <- sort(unique(md[[v]]))
  data.frame(variable = v, level = lv,
             pct_cancer  = sapply(lv, function(l) round(100 * mean(md[[v]][md$cancer == 1] == l), 1)),
             pct_control = sapply(lv, function(l) round(100 * mean(md[[v]][md$cancer == 0] == l), 1)),
             SMD = sapply(lv, function(l) round(smd_bin(md[[v]] == l, md$cancer), 3)), row.names = NULL)
}))
ex_raw$cs_doc <- !is.na(ex_raw$code_status_raw) & ex_raw$code_status_raw != ""
RES$v7_code_doc_by_cancer <- as.data.frame(table(cancer = ex_raw$cancer, documented = ex_raw$cs_doc))
RES$v7_code_doc_by_period <- as.data.frame(table(period = ex_raw$anchor_year_group, documented = ex_raw$cs_doc))
# ---- 15g. 稿子报告的 v7 数字:预期值(主校正集 B)---------------------------------- [V8]
u7 <- function(o, s) v7_uncond[v7_uncond$outcome == o & v7_uncond$set == s, ]
chk_row <- function(prefix, e, r) chk_vec(prefix, e, list(estimate = r$RR, lo = r$RR_lo, hi = r$RR_hi),
                                         list(estimate = r$RD_pts / 100, lo = r$RD_lo / 100, hi = r$RD_hi / 100))
chk_row("v7 无条件 d90 A", c(1.99, 1.90, 2.10, 22.0, 20.3, 23.7), u7("d90", "A"))
chk("v7 E值 点估计 d90 A", u7("d90", "A")$E_point, 3.40)
chk("v7 E值 区间限 d90 A", u7("d90", "A")$E_limit, 3.20)
chk_row("v7 无条件 d90 B", c(1.87, 1.78, 1.97, 19.6, 17.9, 21.4), u7("d90", "B"))
chk("v7 E值 点估计 d90 B", u7("d90", "B")$E_point, 3.15)
chk("v7 E值 区间限 d90 B", u7("d90", "B")$E_limit, 2.95)
chk_row("v7 无条件 d90 C", c(1.86, 1.77, 1.96, 19.4, 17.7, 21.2), u7("d90", "C"))
chk("v7 E值 点估计 d90 C", u7("d90", "C")$E_point, 3.13)
chk("v7 E值 区间限 d90 C", u7("d90", "C")$E_limit, 2.93)
chk_row("v7 无条件 hosp_death A", c(1.62, 1.51, 1.74, 9.2, 7.8, 10.6), u7("hosp_death", "A"))
chk("v7 E值 点估计 hosp_death A", u7("hosp_death", "A")$E_point, 2.62)
chk("v7 E值 区间限 hosp_death A", u7("hosp_death", "A")$E_limit, 2.38)
chk_row("v7 无条件 hosp_death B", c(1.51, 1.40, 1.63, 7.6, 6.1, 9.0), u7("hosp_death", "B"))
chk("v7 E值 点估计 hosp_death B", u7("hosp_death", "B")$E_point, 2.39)
chk("v7 E值 区间限 hosp_death B", u7("hosp_death", "B")$E_limit, 2.16)
chk_row("v7 无条件 hosp_death C", c(1.50, 1.39, 1.61, 7.4, 6.0, 8.9), u7("hosp_death", "C"))
chk("v7 E值 点估计 hosp_death C", u7("hosp_death", "C")$E_point, 2.36)
chk("v7 E值 区间限 hosp_death C", u7("hosp_death", "C")$E_limit, 2.13)
cr_all <- RES$v7_crude_all
chk("v7 粗死亡率 90天% 非肿瘤", 100 * cr_all$d90[cr_all$cancer == 0], 21.6, tol = 0.051)
chk("v7 粗死亡率 住院% 非肿瘤", 100 * cr_all$hosp_death[cr_all$cancer == 0], 14.1, tol = 0.051)
chk("v7 粗死亡率 90天% 肿瘤", 100 * cr_all$d90[cr_all$cancer == 1], 44.4, tol = 0.051)
chk("v7 粗死亡率 住院% 肿瘤", 100 * cr_all$hosp_death[cr_all$cancer == 1], 24.0, tol = 0.051)
chk("v7 逐步缩小样本 n 1", v7_st$n[1], 12636, tol = 0)
chk("v7 逐步缩小样本 RD 1", v7_st$RD_pts[1], 19.6, tol = 0.051)
chk("v7 逐步缩小样本 RR 1", v7_st$RR[1], 1.87)
chk("v7 逐步缩小样本 n 2", v7_st$n[2], 11771, tol = 0)
chk("v7 逐步缩小样本 RD 2", v7_st$RD_pts[2], 19.1, tol = 0.051)
chk("v7 逐步缩小样本 RR 2", v7_st$RR[2], 2.08)
chk("v7 逐步缩小样本 n 3", v7_st$n[3], 8835, tol = 0)
chk("v7 逐步缩小样本 RD 3", v7_st$RD_pts[3], 18.7, tol = 0.051)
chk("v7 逐步缩小样本 RR 3", v7_st$RR[3], 1.65)
chk("v7 逐步缩小样本 n 4", v7_st$n[4], 7970, tol = 0)
chk("v7 逐步缩小样本 RD 4", v7_st$RD_pts[4], 18.3, tol = 0.051)
chk("v7 逐步缩小样本 RR 4", v7_st$RR[4], 1.82)
chk("v7 逐步缩小样本 n 5", v7_st$n[5], 7970, tol = 0)
chk("v7 逐步缩小样本 RD 5", v7_st$RD_pts[5], 17.9, tol = 0.051)
chk("v7 逐步缩小样本 RR 5", v7_st$RR[5], 1.79)
chk("v7 早期出ICU n 2a_left_ICU_<24h 癌症=0", sum(v7_ee$sub == "2a_left_ICU_<24h" & v7_ee$cancer == 0), 2412, tol = 0)
chk("v7 早期出ICU n 2b_left_ICU_24-72h 癌症=0", sum(v7_ee$sub == "2b_left_ICU_24-72h" & v7_ee$cancer == 0), 228, tol = 0)
chk("v7 早期出ICU n 2c_in_ICU_>=72h_no_data 癌症=0", sum(v7_ee$sub == "2c_in_ICU_>=72h_no_data" & v7_ee$cancer == 0), 4, tol = 0)
chk("v7 早期出ICU n 2a_left_ICU_<24h 癌症=1", sum(v7_ee$sub == "2a_left_ICU_<24h" & v7_ee$cancer == 1), 1045, tol = 0)
chk("v7 早期出ICU n 2b_left_ICU_24-72h 癌症=1", sum(v7_ee$sub == "2b_left_ICU_24-72h" & v7_ee$cancer == 1), 110, tol = 0)
chk("v7 早期出ICU n 2c_in_ICU_>=72h_no_data 癌症=1", sum(v7_ee$sub == "2c_in_ICU_>=72h_no_data" & v7_ee$cancer == 1), 2, tol = 0)
ee_fit7 <- RES$v7_ee_fit
chk_row("v7 早期出ICU 2a_left_ICU_<24h", c(3.60, 3.04, 4.25, 20.6, 17.6, 23.7), ee_fit7[ee_fit7$sub == "2a_left_ICU_<24h", ])
chk_row("v7 早期出ICU 2b_left_ICU_24-72h", c(3.38, 2.10, 5.44, 22.9, 13.2, 32.5), ee_fit7[ee_fit7$sub == "2b_left_ICU_24-72h", ])
chk_row("v7 早期出ICU all early ICU exit", c(3.52, 3.01, 4.12, 20.8, 17.9, 23.7), ee_fit7[ee_fit7$sub == "all early ICU exit", ])
ee_sens7 <- RES$v7_ee_sens
chk_row("v7 早期出ICU加回 1", c(1.79, 1.68, 1.92, 17.9, 15.7, 20.0), ee_sens7[1, ])
chk_row("v7 早期出ICU加回 2", c(2.04, 1.92, 2.18, 18.6, 16.9, 20.4), ee_sens7[2, ])
chk_row("v7 早期出ICU加回 3", c(1.80, 1.68, 1.92, 17.9, 15.7, 20.1), ee_sens7[3, ])
chk_row("v7 早期出ICU加回 4", c(2.04, 1.92, 2.18, 18.6, 16.9, 20.4), ee_sens7[4, ])
bal7 <- RES$v7_balance
chk("v7 平衡 adm_type elective/surgical 肿瘤%", bal7$pct_cancer[bal7$variable == "adm_type" & bal7$level == "elective/surgical"], 18.3, tol = 0.051)
chk("v7 平衡 adm_type elective/surgical 对照%", bal7$pct_control[bal7$variable == "adm_type" & bal7$level == "elective/surgical"], 14.1, tol = 0.051)
chk("v7 平衡 adm_type elective/surgical SMD", bal7$SMD[bal7$variable == "adm_type" & bal7$level == "elective/surgical"], 0.114, tol = 0.0006)
chk("v7 平衡 adm_type emergency/urgent 肿瘤%", bal7$pct_cancer[bal7$variable == "adm_type" & bal7$level == "emergency/urgent"], 64.9, tol = 0.051)
chk("v7 平衡 adm_type emergency/urgent 对照%", bal7$pct_control[bal7$variable == "adm_type" & bal7$level == "emergency/urgent"], 71.2, tol = 0.051)
chk("v7 平衡 adm_type emergency/urgent SMD", bal7$SMD[bal7$variable == "adm_type" & bal7$level == "emergency/urgent"], -0.136, tol = 0.0006)
chk("v7 平衡 adm_type observation/other 肿瘤%", bal7$pct_cancer[bal7$variable == "adm_type" & bal7$level == "observation/other"], 16.8, tol = 0.051)
chk("v7 平衡 adm_type observation/other 对照%", bal7$pct_control[bal7$variable == "adm_type" & bal7$level == "observation/other"], 14.7, tol = 0.051)
chk("v7 平衡 adm_type observation/other SMD", bal7$SMD[bal7$variable == "adm_type" & bal7$level == "observation/other"], 0.058, tol = 0.0006)
chk("v7 平衡 unit Cardiac 肿瘤%", bal7$pct_cancer[bal7$variable == "unit" & bal7$level == "Cardiac"], 11.6, tol = 0.051)
chk("v7 平衡 unit Cardiac 对照%", bal7$pct_control[bal7$variable == "unit" & bal7$level == "Cardiac"], 32.4, tol = 0.051)
chk("v7 平衡 unit Cardiac SMD", bal7$SMD[bal7$variable == "unit" & bal7$level == "Cardiac"], -0.521, tol = 0.0006)
chk("v7 平衡 unit MICU 肿瘤%", bal7$pct_cancer[bal7$variable == "unit" & bal7$level == "MICU"], 16.8, tol = 0.051)
chk("v7 平衡 unit MICU 对照%", bal7$pct_control[bal7$variable == "unit" & bal7$level == "MICU"], 19.6, tol = 0.051)
chk("v7 平衡 unit MICU SMD", bal7$SMD[bal7$variable == "unit" & bal7$level == "MICU"], -0.073, tol = 0.0006)
chk("v7 平衡 unit MICU/SICU 肿瘤%", bal7$pct_cancer[bal7$variable == "unit" & bal7$level == "MICU/SICU"], 36.0, tol = 0.051)
chk("v7 平衡 unit MICU/SICU 对照%", bal7$pct_control[bal7$variable == "unit" & bal7$level == "MICU/SICU"], 9.6, tol = 0.051)
chk("v7 平衡 unit MICU/SICU SMD", bal7$SMD[bal7$variable == "unit" & bal7$level == "MICU/SICU"], 0.662, tol = 0.0006)
chk("v7 平衡 unit Neuro 肿瘤%", bal7$pct_cancer[bal7$variable == "unit" & bal7$level == "Neuro"], 7.6, tol = 0.051)
chk("v7 平衡 unit Neuro 对照%", bal7$pct_control[bal7$variable == "unit" & bal7$level == "Neuro"], 11.7, tol = 0.051)
chk("v7 平衡 unit Neuro SMD", bal7$SMD[bal7$variable == "unit" & bal7$level == "Neuro"], -0.139, tol = 0.0006)
chk("v7 平衡 unit Surgical/trauma/other 肿瘤%", bal7$pct_cancer[bal7$variable == "unit" & bal7$level == "Surgical/trauma/other"], 28.1, tol = 0.051)
chk("v7 平衡 unit Surgical/trauma/other 对照%", bal7$pct_control[bal7$variable == "unit" & bal7$level == "Surgical/trauma/other"], 26.7, tol = 0.051)
chk("v7 平衡 unit Surgical/trauma/other SMD", bal7$SMD[bal7$variable == "unit" & bal7$level == "Surgical/trauma/other"], 0.031, tol = 0.0006)
chk("v7 平衡 svc medical/other 肿瘤%", bal7$pct_cancer[bal7$variable == "svc" & bal7$level == "medical/other"], 65.0, tol = 0.051)
chk("v7 平衡 svc medical/other 对照%", bal7$pct_control[bal7$variable == "svc" & bal7$level == "medical/other"], 58.1, tol = 0.051)
chk("v7 平衡 svc medical/other SMD", bal7$SMD[bal7$variable == "svc" & bal7$level == "medical/other"], 0.142, tol = 0.0006)
chk("v7 平衡 svc surgical 肿瘤%", bal7$pct_cancer[bal7$variable == "svc" & bal7$level == "surgical"], 35.0, tol = 0.051)
chk("v7 平衡 svc surgical 对照%", bal7$pct_control[bal7$variable == "svc" & bal7$level == "surgical"], 41.9, tol = 0.051)
chk("v7 平衡 svc surgical SMD", bal7$SMD[bal7$variable == "svc" & bal7$level == "surgical"], -0.142, tol = 0.0006)
chk("v7 平衡 period 2008 - 2010 肿瘤%", bal7$pct_cancer[bal7$variable == "period" & bal7$level == "2008 - 2010"], 27.8, tol = 0.051)
chk("v7 平衡 period 2008 - 2010 对照%", bal7$pct_control[bal7$variable == "period" & bal7$level == "2008 - 2010"], 25.1, tol = 0.051)
chk("v7 平衡 period 2008 - 2010 SMD", bal7$SMD[bal7$variable == "period" & bal7$level == "2008 - 2010"], 0.061, tol = 0.0006)
chk("v7 平衡 period 2011 - 2013 肿瘤%", bal7$pct_cancer[bal7$variable == "period" & bal7$level == "2011 - 2013"], 20.2, tol = 0.051)
chk("v7 平衡 period 2011 - 2013 对照%", bal7$pct_control[bal7$variable == "period" & bal7$level == "2011 - 2013"], 19.5, tol = 0.051)
chk("v7 平衡 period 2011 - 2013 SMD", bal7$SMD[bal7$variable == "period" & bal7$level == "2011 - 2013"], 0.018, tol = 0.0006)
chk("v7 平衡 period 2014 - 2016 肿瘤%", bal7$pct_cancer[bal7$variable == "period" & bal7$level == "2014 - 2016"], 20.0, tol = 0.051)
chk("v7 平衡 period 2014 - 2016 对照%", bal7$pct_control[bal7$variable == "period" & bal7$level == "2014 - 2016"], 20.7, tol = 0.051)
chk("v7 平衡 period 2014 - 2016 SMD", bal7$SMD[bal7$variable == "period" & bal7$level == "2014 - 2016"], -0.018, tol = 0.0006)
chk("v7 平衡 period 2017 - 2019 肿瘤%", bal7$pct_cancer[bal7$variable == "period" & bal7$level == "2017 - 2019"], 18.5, tol = 0.051)
chk("v7 平衡 period 2017 - 2019 对照%", bal7$pct_control[bal7$variable == "period" & bal7$level == "2017 - 2019"], 19.4, tol = 0.051)
chk("v7 平衡 period 2017 - 2019 SMD", bal7$SMD[bal7$variable == "period" & bal7$level == "2017 - 2019"], -0.025, tol = 0.0006)
chk("v7 平衡 period 2020 - 2022 肿瘤%", bal7$pct_cancer[bal7$variable == "period" & bal7$level == "2020 - 2022"], 13.6, tol = 0.051)
chk("v7 平衡 period 2020 - 2022 对照%", bal7$pct_control[bal7$variable == "period" & bal7$level == "2020 - 2022"], 15.3, tol = 0.051)
chk("v7 平衡 period 2020 - 2022 SMD", bal7$SMD[bal7$variable == "period" & bal7$level == "2020 - 2022"], -0.050, tol = 0.0006)
chk("v7 平衡 code comfort measures 肿瘤%", bal7$pct_cancer[bal7$variable == "code" & bal7$level == "comfort measures"], 0.2, tol = 0.051)
chk("v7 平衡 code comfort measures 对照%", bal7$pct_control[bal7$variable == "code" & bal7$level == "comfort measures"], 0.2, tol = 0.051)
chk("v7 平衡 code comfort measures SMD", bal7$SMD[bal7$variable == "code" & bal7$level == "comfort measures"], -0.003, tol = 0.0006)
chk("v7 平衡 code DNR/DNI 肿瘤%", bal7$pct_cancer[bal7$variable == "code" & bal7$level == "DNR/DNI"], 3.5, tol = 0.051)
chk("v7 平衡 code DNR/DNI 对照%", bal7$pct_control[bal7$variable == "code" & bal7$level == "DNR/DNI"], 2.6, tol = 0.051)
chk("v7 平衡 code DNR/DNI SMD", bal7$SMD[bal7$variable == "code" & bal7$level == "DNR/DNI"], 0.053, tol = 0.0006)
chk("v7 平衡 code full code 肿瘤%", bal7$pct_cancer[bal7$variable == "code" & bal7$level == "full code"], 27.1, tol = 0.051)
chk("v7 平衡 code full code 对照%", bal7$pct_control[bal7$variable == "code" & bal7$level == "full code"], 36.0, tol = 0.051)
chk("v7 平衡 code full code SMD", bal7$SMD[bal7$variable == "code" & bal7$level == "full code"], -0.193, tol = 0.0006)
chk("v7 平衡 code not documented 肿瘤%", bal7$pct_cancer[bal7$variable == "code" & bal7$level == "not documented"], 69.3, tol = 0.051)
chk("v7 平衡 code not documented 对照%", bal7$pct_control[bal7$variable == "code" & bal7$level == "not documented"], 61.3, tol = 0.051)
chk("v7 平衡 code not documented SMD", bal7$SMD[bal7$variable == "code" & bal7$level == "not documented"], 0.169, tol = 0.0006)
chk("v7 code status 记录 癌症=0 有记录=False", sum(ex_raw$cancer == 0 & ex_raw$cs_doc == FALSE), 11766, tol = 0)
chk("v7 code status 记录 癌症=1 有记录=False", sum(ex_raw$cancer == 1 & ex_raw$cs_doc == FALSE), 2919, tol = 0)
chk("v7 code status 记录 癌症=0 有记录=True", sum(ex_raw$cancer == 0 & ex_raw$cs_doc == TRUE), 7569, tol = 0)
chk("v7 code status 记录 癌症=1 有记录=True", sum(ex_raw$cancer == 1 & ex_raw$cs_doc == TRUE), 1294, tol = 0)
chk("v7 code status 记录 2008 - 2010 有记录=False", sum(ex_raw$anchor_year_group == "2008 - 2010" & ex_raw$cs_doc == FALSE), 3096, tol = 0)
chk("v7 code status 记录 2011 - 2013 有记录=False", sum(ex_raw$anchor_year_group == "2011 - 2013" & ex_raw$cs_doc == FALSE), 2527, tol = 0)
chk("v7 code status 记录 2014 - 2016 有记录=False", sum(ex_raw$anchor_year_group == "2014 - 2016" & ex_raw$cs_doc == FALSE), 3086, tol = 0)
chk("v7 code status 记录 2017 - 2019 有记录=False", sum(ex_raw$anchor_year_group == "2017 - 2019" & ex_raw$cs_doc == FALSE), 3156, tol = 0)
chk("v7 code status 记录 2020 - 2022 有记录=False", sum(ex_raw$anchor_year_group == "2020 - 2022" & ex_raw$cs_doc == FALSE), 2820, tol = 0)
chk("v7 code status 记录 2008 - 2010 有记录=True", sum(ex_raw$anchor_year_group == "2008 - 2010" & ex_raw$cs_doc == TRUE), 2987, tol = 0)
chk("v7 code status 记录 2011 - 2013 有记录=True", sum(ex_raw$anchor_year_group == "2011 - 2013" & ex_raw$cs_doc == TRUE), 2090, tol = 0)
chk("v7 code status 记录 2014 - 2016 有记录=True", sum(ex_raw$anchor_year_group == "2014 - 2016" & ex_raw$cs_doc == TRUE), 1621, tol = 0)
chk("v7 code status 记录 2017 - 2019 有记录=True", sum(ex_raw$anchor_year_group == "2017 - 2019" & ex_raw$cs_doc == TRUE), 1335, tol = 0)
chk("v7 code status 记录 2020 - 2022 有记录=True", sum(ex_raw$anchor_year_group == "2020 - 2022" & ex_raw$cs_doc == TRUE), 830, tol = 0)

cat("\n[v7] 第 15 节完成:结果表以 v7_ 开头写入 reproduce_output/\n")

# ---- 14. 输出 --------------------------------------------------------------------
write.csv(CHK, file.path(OUT_DIR, "checks.csv"), row.names = FALSE)
for (nm in names(RES)) {
  x <- RES[[nm]]
  if (is.list(x) && !is.data.frame(x)) {
    for (k in names(x)) write.csv(x[[k]], file.path(OUT_DIR, paste0(nm, "_", k, ".csv")))
  } else write.csv(x, file.path(OUT_DIR, paste0(nm, ".csv")))
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "session_info.txt"))

cat("\n==== 汇总 ====\n")
print(table(CHK$status))
if (any(CHK$status == "DIFF")) {
  cat("\n---- 与稿子不一致的项 ----\n")
  print(CHK[CHK$status == "DIFF", ], row.names = FALSE)
}
cat("\n输出目录:", OUT_DIR, "\n")
