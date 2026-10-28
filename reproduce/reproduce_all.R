# =============================================================================
# reproduce_all.R
# Trajectories of acute kidney injury and 90-day mortality in critically ill
# patients with and without cancer (MIMIC-IV v3.1)
#
# 复现稿子中的 Table 1–7、补充图 S1 的数字、AUC,以及 RRT / 安宁疗护的核查。
#
# 每个检查项都带状态标记:
#   [R]  2026-10-23 已在最终版数据上完全复现
#   [R*] 已复现,但与稿子有 0.1–0.4 个百分点的细微出入(原因未查明,见 note 列)
#   [N]  本脚本新增,首次运行,结果未经核对
#
# 运行前:
#   1. 在仍保留分析对象的会话里运行 00_export_inputs.R(生成 inputs_final.rds)
#   2. 修改下面的 MIMIC_DIR(第 9、10 节需要 admissions / procedureevents)
# 运行: source("reproduce_all.R", echo = FALSE)
# 输出: reproduce_output/ 下的 checks.csv(逐项对照稿子)、各表 csv、session_info.txt
# =============================================================================

# ---- 0. 配置 -----------------------------------------------------------------
if (!exists("PROJ_DIR")) PROJ_DIR <- getwd()
# 下面四个变量,如果在运行 source() 之前已经在控制台里设好,就沿用控制台里的值:
#   ADM_PATH <- file.choose()   # 弹出对话框,选 admissions.csv.gz
#   PE_PATH  <- file.choose()   # 选 procedureevents.csv.gz
#   DI_PATH  <- file.choose()   # 选 d_items.csv.gz(可选;没有就设为 NA_character_)
# 或者设 MIMIC_DIR(会在该文件夹内递归查找;不要设成主目录,否则查找会非常慢)
if (!exists("MIMIC_DIR")) MIMIC_DIR <- NA_character_
if (!exists("ADM_PATH"))  ADM_PATH  <- NA_character_   # admissions.csv.gz(取 discharge_location)
if (!exists("PE_PATH"))   PE_PATH   <- NA_character_   # procedureevents.csv.gz
if (!exists("DI_PATH"))   DI_PATH   <- NA_character_   # d_items.csv.gz(可选)
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

cov_names <- c("age", "male", "sofa", "sepsis", "chf", "ckd", "copd", "dm", "liver")
covs <- paste(cov_names, collapse = " + ")
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
num  <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else as.numeric(x)

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

# 肿瘤 vs 非肿瘤的校正 OR(按匹配集聚类的稳健标准误)
or_row <- function(dat, outcome, extra = "") {
  f   <- glm(as.formula(paste(outcome, "~ cancer +", covs, extra)),
             data = dat, family = binomial)
  vc  <- vcovCL(f, cluster = dat$subclass)
  est <- coef(f)["cancer"]; se <- sqrt(vc["cancer", "cancer"])
  c(n = nrow(dat), OR = unname(exp(est)),
    lo = unname(exp(est - 1.96 * se)), hi = unname(exp(est + 1.96 * se)))
}
run_or <- function(dat, y) {
  parts <- lapply(split(dat, dat$pheno6), or_row, outcome = y)
  rbind(do.call(rbind, parts), All = or_row(dat, y, "+ factor(pheno6)"))
}

# 肾脏结局(口径按稿子的文字描述写成;见 [R*] 说明)
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

# Table 6:各组相对"峰值 2 期 / 一过性"的校正 OR
or_grp <- function(dat, y) {
  d <- dat[!is.na(dat[[y]]), ]
  f <- glm(as.formula(paste(y, "~ factor(g) + cancer +", covs)), data = d, family = binomial)
  ct <- coeftest(f, vcov = vcovCL(f, cluster = d$subclass))
  idx <- grep("factor\\(g\\)", rownames(ct))
  est <- ct[idx, 1]; se <- ct[idx, 2]
  out <- cbind(OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se))
  rownames(out) <- sub("factor\\(g\\)", "", rownames(out))
  out
}
# Table 7:肿瘤 vs 非肿瘤(含事件数)
or_c <- function(dat, y) {
  d <- dat[!is.na(dat[[y]]), ]
  r <- or_row(d, y, "+ factor(pheno6)")
  c(r, n_c = sum(d$cancer == 1), ev_c = sum(d[[y]][d$cancer == 1]),
       n_k = sum(d$cancer == 0), ev_k = sum(d[[y]][d$cancer == 0]))
}
# 各轨迹相对一过性的校正 OR(队列内)
within_or <- function(dat, y) {
  d <- dat
  d$tr <- factor(d$pheno6, levels = c("3_transient", "4_persistent_stage2",
                                      "5_worsening", "6_severe_from_start"))
  f  <- glm(as.formula(paste(y, "~ tr +", covs)), data = d, family = binomial)
  ix <- grep("^tr", names(coef(f)))
  ci <- confint.default(f)[ix, , drop = FALSE]
  out <- cbind(OR = exp(coef(f)[ix]), lo = exp(ci[, 1]), hi = exp(ci[, 2]))
  rownames(out) <- sub("^tr", "", rownames(out))
  out
}
wald_p <- function(f0, f1, cl) {
  w <- waldtest(f0, f1, vcov = function(x) vcovCL(x, cluster = cl), test = "Chisq")
  c(chisq = w$Chisq[2], df = w$Df[2], p = w[["Pr(>Chisq)"]][2])
}

# 读取 gz 压缩 csv 的指定列(只读需要的列,省内存)
read_gz_cols <- function(path, cols) {
  hdr <- gsub('"', "", strsplit(readLines(gzfile(path), n = 1), ",")[[1]])
  cc  <- ifelse(hdr %in% cols, "character", "NULL"); names(cc) <- hdr
  read.csv(gzfile(path), colClasses = cc, na.strings = c("", "NA"))
}

# 逐项对照稿子
CHK <- NULL
chk <- function(item, observed, expected, tol = 0.0051, note = "") {
  ok <- length(observed) == 1 && is.finite(observed) && abs(observed - expected) <= tol
  CHK <<- rbind(CHK, data.frame(item = item, expected = expected,
                                observed = round(as.numeric(observed), 4),
                                status = if (ok) "OK" else "DIFF", note = note,
                                stringsAsFactors = FALSE))
  invisible(ok)
}
mk <- function(rn, OR, lo, hi) matrix(c(OR, lo, hi), ncol = 3,
                                      dimnames = list(rn, c("OR", "lo", "hi")))
chk_or <- function(prefix, res, exp, note = "") {
  for (r in rownames(exp)) for (k in c("OR", "lo", "hi"))
    chk(paste(prefix, r, k), res[r, k], exp[r, k], note = note)
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

# ---- 4. Table 2:匹配后平衡(SMD)------------------------------------------- [N]
smd_fun <- function(x, g) {
  x1 <- x[g == 1]; x0 <- x[g == 0]
  (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
}
t2 <- data.frame(variable = cov_names,
                 smd = sapply(cov_names, function(v) smd_fun(num(md[[v]]), md$cancer)))
exp_smd <- c(age = .015, male = .004, sofa = .044, sepsis = .001, chf = .040,
             ckd = .038, copd = .033, dm = .026, liver = .038)
for (v in cov_names) chk(paste("Table 2 SMD", v), abs(t2$smd[t2$variable == v]),
                         exp_smd[[v]], tol = 0.0006)
RES$table2 <- t2

# ---- 5. 主分析样本;Table 4 与交互检验 ---------------------------------------- [R]
md3 <- md |>
  select(stay_id, cancer, subclass, all_of(cov_names), hosp_death, d90) |>
  left_join(traj_all |> select(stay_id, pheno6), by = "stay_id") |>
  filter(!pheno6 %in% c("1_early_death", "2_early_icu_exit"))
chk("主样本 n", nrow(md3), 7970, tol = 0)
chk("主样本 肿瘤 n", sum(md3$cancer == 1), 2656, tol = 0)
chk("主样本 非肿瘤 n", sum(md3$cancer == 0), 5314, tol = 0)

rn4 <- c("3_transient", "4_persistent_stage2", "5_worsening", "6_severe_from_start", "All")
t4_90 <- run_or(md3, "d90")
t4_h  <- run_or(md3, "hosp_death")
chk_or("Table 4 d90", t4_90, mk(rn4, c(3.10, 3.23, 2.14, 2.35, 2.84),
                                     c(2.58, 2.68, 1.50, 1.86, 2.54),
                                     c(3.74, 3.90, 3.06, 2.96, 3.16)))
chk_or("Table 4 hosp", t4_h, mk(rn4, c(2.02, 2.05, 1.19, 1.40, 1.70),
                                     c(1.57, 1.63, 0.82, 1.09, 1.49),
                                     c(2.62, 2.58, 1.72, 1.79, 1.95)))
RES$table4_d90 <- as.data.frame(t4_90); RES$table4_hosp <- as.data.frame(t4_h)

fm <- function(y, rhs, fn = glm) {
  if (identical(fn, glm)) glm(as.formula(paste(y, "~", rhs)), data = md3, family = binomial)
  else lm(as.formula(paste(y, "~", rhs)), data = md3)
}
w90  <- wald_p(fm("d90", paste("cancer + factor(pheno6) +", covs)),
               fm("d90", paste("cancer * factor(pheno6) +", covs)), md3$subclass)
whos <- wald_p(fm("hosp_death", paste("cancer + factor(pheno6) +", covs)),
               fm("hosp_death", paste("cancer * factor(pheno6) +", covs)), md3$subclass)
wadd <- wald_p(fm("d90", paste("cancer + factor(pheno6) +", covs), lm),
               fm("d90", paste("cancer * factor(pheno6) +", covs), lm), md3$subclass)
chk("交互 乘法 d90 chisq", w90["chisq"], 9.7449, tol = 0.001)
chk("交互 乘法 d90 P", w90["p"], 0.021, tol = 0.0006)
chk("交互 乘法 hosp chisq", whos["chisq"], 11.257, tol = 0.001)
chk("交互 乘法 hosp P", whos["p"], 0.010, tol = 0.0006)
chk("交互 加法 d90 chisq", wadd["chisq"], 4.1442, tol = 0.001)
chk("交互 加法 d90 P", wadd["p"], 0.25, tol = 0.0051)

# 各队列内、相对一过性的 OR(landmark 分析,稿子"结果"第二节)--------------------- [N]
wc <- md3 |> filter(cancer == 1); wk <- md3 |> filter(cancer == 0)
r3 <- c("4_persistent_stage2", "5_worsening", "6_severe_from_start")
chk_or("队列内 肿瘤 d90", within_or(wc, "d90"),
       mk(r3, c(1.47, 2.18, 1.89), c(1.21, 1.62, 1.51), c(1.79, 2.93, 2.37)))
chk_or("队列内 肿瘤 hosp", within_or(wc, "hosp_death"),
       mk(r3, c(1.72, 2.76, 2.01), c(1.33, 1.96, 1.52), c(2.22, 3.89, 2.66)))
chk_or("队列内 非肿瘤 d90", within_or(wk, "d90"),
       mk(r3, c(1.35, 3.04, 2.23), c(1.13, 2.35, 1.84), c(1.61, 3.92, 2.71)))

# 补充图 S1:加回"早期出 ICU"类 ---------------------------------------------- [N]
md_alive <- md |>
  select(stay_id, cancer, subclass, all_of(cov_names), hosp_death, d90) |>
  left_join(traj_all |> select(stay_id, pheno6), by = "stay_id") |>
  filter(pheno6 != "1_early_death")
s1 <- md_alive |> filter(pheno6 == "2_early_icu_exit")
chk("S1 早期出ICU n", nrow(s1), 3801, tol = 0)
chk("S1 早期出ICU d90 OR", or_row(s1, "d90")["OR"], 5.86)
chk("S1 早期出ICU d90 lo", or_row(s1, "d90")["lo"], 4.82)
chk("S1 早期出ICU d90 hi", or_row(s1, "d90")["hi"], 7.12)
chk("S1 早期出ICU hosp OR", or_row(s1, "hosp_death")["OR"], 5.62)
chk("S1 早期出ICU hosp lo", or_row(s1, "hosp_death")["lo"], 3.70)
chk("S1 早期出ICU hosp hi", or_row(s1, "hosp_death")["hi"], 8.54)
chk("S1 加回后总体 d90 OR", or_row(md_alive, "d90", "+ factor(pheno6)")["OR"], 3.38)

# ---- 6. Table 5:按肿瘤类型 --------------------------------------------------- [R]
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
r5 <- c("Hematologic", "Non-metastatic solid", "Metastatic solid")
t5_h  <- do.call(rbind, lapply(split(md4, md4$ctype), or_row, outcome = "hosp_death"))
t5_90 <- do.call(rbind, lapply(split(md4, md4$ctype), or_row, outcome = "d90"))
chk_or("Table 5 hosp", t5_h, mk(r5, c(2.31, 0.88, 2.65), c(1.73, 0.70, 2.17), c(3.08, 1.11, 3.23)))
chk_or("Table 5 d90",  t5_90, mk(r5, c(2.89, 1.35, 5.67), c(2.25, 1.12, 4.78), c(3.72, 1.62, 6.74)))
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
het <- function(y) {
  f0 <- glm(as.formula(paste(y, "~ cancer + factor(ctype) +", covs)), data = md4, family = binomial)
  f1 <- glm(as.formula(paste(y, "~ cancer * factor(ctype) +", covs)), data = md4, family = binomial)
  wald_p(f0, f1, md4$subclass)
}
chk("Table 5 异质性 hosp chisq", het("hosp_death")["chisq"], 52.437, tol = 0.01)
chk("Table 5 异质性 d90 chisq",  het("d90")["chisq"], 130.66, tol = 0.01)
RES$table5_hosp <- as.data.frame(t5_h); RES$table5_d90 <- as.data.frame(t5_90)

# ---- 7. Table 6:峰值分期 × 轨迹;肾脏结局 ---------------------------------- [R]/[R*]
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
G <- c("A_2_transient", "B_3_resolved", "C_2_persistent", "D_3_worsening", "E_3_severe")

s6 <- grp2 |> group_by(g) |>
  summarise(n = n(), d90 = 100 * mean(d90),
            rrt_late = 100 * mean(y_rrt_late, na.rm = TRUE),
            nonrec = 100 * mean(y_nonrec, na.rm = TRUE),
            make = 100 * mean(y_make), .groups = "drop") |> as.data.frame()
rownames(s6) <- s6$g
nt <- "已知小差异:口径按稿子文字重写,原因未查明"
for (i in seq_along(G)) {
  chk(paste("Table 6 n", G[i]), s6[G[i], "n"], c(2907, 276, 2626, 620, 1541)[i], tol = 0)
  chk(paste("Table 6 d90%", G[i]), s6[G[i], "d90"], c(18.9, 27.9, 27.0, 45.0, 42.4)[i], tol = 0.051)
  chk(paste("Table 6 72h后RRT%", G[i]), s6[G[i], "rrt_late"], c(0.6, 2.3, 1.6, 8.4, 6.7)[i], tol = 0.051)
  chk(paste("Table 6 非恢复%", G[i]), s6[G[i], "nonrec"], c(6.8, 9.2, 11.2, 24.0, 29.4)[i], tol = 0.051, note = nt)
  chk(paste("Table 6 MAKE%", G[i]), s6[G[i], "make"], c(10.1, 17.0, 19.1, 43.9, 44.1)[i], tol = 0.051, note = nt)
}
o_d90 <- or_grp(grp2, "d90"); o_late <- or_grp(grp2, "y_rrt_late")
o_any <- or_grp(grp2, "y_rrt_any"); o_nr <- or_grp(grp2, "y_nonrec"); o_mk <- or_grp(grp2, "y_make")
r6 <- c("B_3_resolved", "C_2_persistent", "D_3_worsening", "E_3_severe")
chk_or("Table 6 d90 OR", o_d90, mk(r6, c(1.48, 1.47, 2.76, 2.17), c(1.11, 1.28, 2.26, 1.87), c(1.98, 1.68, 3.37, 2.52)))
chk_or("Table 6 72h后RRT OR", o_late, mk(r6, c(3.30, 2.66, 10.86, 6.75), c(1.29, 1.50, 6.09, 3.95), c(8.42, 4.71, 19.37, 11.54)))
for (i in seq_along(r6)) {
  chk(paste("Table 6 任何RRT OR", r6[i]), o_any[r6[i], "OR"], c(6.41, 2.36, 16.20, 15.27)[i])
  chk(paste("Table 6 非恢复 OR", r6[i]), o_nr[r6[i], "OR"], c(1.35, 1.70, 3.69, 3.99)[i], note = nt)
  chk(paste("Table 6 MAKE OR", r6[i]), o_mk[r6[i], "OR"], c(1.62, 2.01, 5.43, 4.42)[i], note = nt)
}
RES$table6_summary <- s6
RES$table6_or <- list(d90 = o_d90, rrt_late = o_late, rrt_any = o_any, nonrec = o_nr, make = o_mk) |>
  lapply(function(m) as.data.frame(round(m, 2)))

# 峰值分期之外,轨迹是否带来额外信息(Wald 检验)------------------------------------ [N]
grp$peak <- ifelse(grp$g %in% c("A_2_transient", "C_2_persistent"), "2", "3")
# 模型 1:峰值分期;模型 2:再加入峰值内的轨迹对比(持续性 vs 一过性;
# 恶化、起病即重度 vs 已缓解),共 3 个自由度。两个模型嵌套,可用 Wald 检验比较。
grp$peak3 <- as.numeric(grp$peak == "3")
grp$pC <- as.numeric(grp$g == "C_2_persistent")
grp$pD <- as.numeric(grp$g == "D_3_worsening")
grp$pE <- as.numeric(grp$g == "E_3_severe")
wpk <- wald_p(glm(as.formula(paste("d90 ~ cancer + peak3 +", covs)), data = grp, family = binomial),
              glm(as.formula(paste("d90 ~ cancer + peak3 + pC + pD + pE +", covs)),
                  data = grp, family = binomial),
              grp$subclass)
chk("轨迹加入峰值分期 Wald chisq", wpk["chisq"], 45.6, tol = 0.06)
chk("轨迹加入峰值分期 df", wpk["df"], 3, tol = 0)

# ---- 8. Table 7:肿瘤 vs 非肿瘤的肾脏结局 ------------------------------------ [R]/[R*]
kd <- md3 |> left_join(kid, by = "stay_id")
t7 <- rbind(any_rrt = or_c(kd, "y_rrt_any"), rrt_late = or_c(kd, "y_rrt_late"),
            nonrec = or_c(kd, "y_nonrec"), make = or_c(kd, "y_make"))
ex7 <- list(any_rrt  = c(.83, .66, 1.05, 2656, 149, 5314, 331),
            rrt_late = c(.80, .57, 1.11, 2569, 62, 5124, 141),
            nonrec   = c(1.15, .98, 1.36, 1895, 274, 4161, 526),
            make     = c(1.62, 1.44, 1.83, 2656, 740, 5314, 1053))
for (k in names(ex7)) {
  nn <- if (k %in% c("nonrec", "make")) nt else ""
  chk(paste("Table 7", k, "OR"), t7[k, "OR"], ex7[[k]][1], note = nn)
  chk(paste("Table 7", k, "lo"), t7[k, "lo"], ex7[[k]][2], note = nn)
  chk(paste("Table 7", k, "hi"), t7[k, "hi"], ex7[[k]][3], note = nn)
  chk(paste("Table 7", k, "n 肿瘤"), t7[k, "n_c"], ex7[[k]][4], tol = 0, note = nn)
  chk(paste("Table 7", k, "事件 肿瘤"), t7[k, "ev_c"], ex7[[k]][5], tol = 0, note = nn)
  chk(paste("Table 7", k, "n 非肿瘤"), t7[k, "n_k"], ex7[[k]][6], tol = 0, note = nn)
  chk(paste("Table 7", k, "事件 非肿瘤"), t7[k, "ev_k"], ex7[[k]][7], tol = 0, note = nn)
}
RES$table7 <- as.data.frame(round(t7, 3))

# ---- 9. AUC(五折交叉验证,折按匹配集分配)------------------------------------ [R]
set.seed(1)
ss <- unique(grp$subclass)
fold_of <- setNames(sample(rep(1:5, length.out = length(ss))), ss)
grp$fold <- fold_of[as.character(grp$subclass)]
forms <- list(m0 = paste("d90 ~ cancer +", covs),
              m1 = paste("d90 ~ cancer + factor(peak) +", covs),
              m2 = paste("d90 ~ cancer + factor(g) +", covs))
cvpred <- function(f) {
  p <- numeric(nrow(grp))
  for (k in 1:5) {
    fit <- glm(as.formula(f), data = grp[grp$fold != k, ], family = binomial)
    p[grp$fold == k] <- predict(fit, newdata = grp[grp$fold == k, ], type = "response")
  }
  p
}
rocs <- lapply(lapply(forms, cvpred), function(p) roc(grp$d90, p, quiet = TRUE))
aucs <- sapply(rocs, function(r) as.numeric(auc(r)))
dl   <- roc.test(rocs$m1, rocs$m2, method = "delong")
chk("AUC 无峰值分期/轨迹", aucs["m0"], 0.728, tol = 0.0006)
chk("AUC 加峰值分期",     aucs["m1"], 0.737, tol = 0.0006)
chk("AUC 加峰值分期+轨迹", aucs["m2"], 0.741, tol = 0.0006)
chk("AUC DeLong P", dl$p.value, 0.023, tol = 0.0006)
RES$auc <- data.frame(model = names(aucs), auc = round(aucs, 4))

# ---- 10. Table 1 要点(肿瘤队列,各轨迹类的事件数)----------------------------- [N]
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
cat("\n---- rrt_h 取整方式(rrt_time 与 t0 之差,共", nrow(rt), "例)----\n",
    "  round  :", round(mean(round(h_raw) == rt$rrt_h), 4), "\n",
    "  floor  :", round(mean(floor(h_raw) == rt$rrt_h), 4), "\n",
    "  ceiling:", round(mean(ceiling(h_raw) == rt$rrt_h), 4), "\n", sep = "")

# ---- 12. 出院去向(安宁疗护)敏感性分析 ---------------------------------------- [N]
find_mimic <- function(explicit, stem) {
  if (!is.na(explicit) && file.exists(explicit)) return(explicit)
  if (!dir.exists(MIMIC_DIR)) return(NA_character_)
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
  s7 <- rbind(d90 = or_c(kdn, "d90"), make = or_c(kdn, "y_make"), nonrec = or_c(kdn, "y_nonrec"))
  exs <- list(d90 = c(2.55, 2.27, 2.86, 2416, 900, 5182, 1039),
              make = c(1.69, 1.49, 1.91, 2416, 692, 5182, 1035),
              nonrec = c(1.07, 0.90, 1.27, 1721, 227, 4061, 501))
  for (k in names(exs)) {
    nn <- if (k == "d90") "" else nt
    chk(paste("排除安宁疗护", k, "OR"), s7[k, "OR"], exs[[k]][1], note = nn)
    chk(paste("排除安宁疗护", k, "lo"), s7[k, "lo"], exs[[k]][2], note = nn)
    chk(paste("排除安宁疗护", k, "hi"), s7[k, "hi"], exs[[k]][3], note = nn)
    chk(paste("排除安宁疗护", k, "事件 肿瘤"), s7[k, "ev_c"], exs[[k]][5], tol = 0, note = nn)
    chk(paste("排除安宁疗护", k, "事件 非肿瘤"), s7[k, "ev_k"], exs[[k]][7], tol = 0, note = nn)
  }
  RES$hospice_survivors <- as.data.frame(surv); RES$table7_no_hospice <- as.data.frame(round(s7, 3))
} else message("【跳过第 12 节】未找到 admissions.csv(.gz),安宁疗护敏感性分析未运行")

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
