# ============================================================
# 08  论文用：Table 1（各轨迹的基线特征）+ 图1（每日分期）+ 图2（生存曲线）
# 【更新版】轨迹分类只用第0~2天（0~71小时），与72小时界标完全一致。
# 前提：R 会话里已有 baseline、outc；daily、dat 在会话里，或在工作目录有 .rds
# 产出（保存在当前工作目录，用 getwd() 查看）：
#     table1_by_trajectory.csv   汇总表（只有汇总数字）
#     fig1_daily_stage.png       每一类轨迹每天的分期构成
#     fig2_km_landmark.png       72小时界标后的生存曲线
# 注意：图1仍显示第0~6天的每日分期，仅用于描述；第3~6天没有参与分类。
# The console prints summary numbers.
# 只用 R 自带的包（dplyr 你已经装过），不需要再安装别的。
# ============================================================

library(dplyr)
library(survival)   # R 自带

safe_max <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }
to_time  <- function(x) if (inherits(x, "POSIXt")) x else as.POSIXct(as.character(x), tz = "UTC")

if (!exists("daily")) daily <- readRDS("daily_stage.rds")
if (!exists("dat"))   dat   <- readRDS("traj_rule_based.rds")

# ---------- 1. 重建完整数据（4213人），只用第0~2天（0~71小时）划分轨迹 ----------
# 修改说明：day = floor(hours_from_t0 / 24)，第3天是72~95小时，已在72小时界标之后。
#           为避免用到界标之后的信息，分类只用第1~2天（加第0天），
#           因此下面三行由 1:3 改为 1:2。变量名 peak13 保留以免改动后面的代码，
#           含义现在是“第1~2天最高分期”。
d72 <- daily %>%
  group_by(stay_id) %>%
  summarise(
    day0_s  = safe_max(max_stage[day == 0]),
    obs     = sum(day %in% 1:2 & !is.na(max_stage)),
    persist = sum(day %in% 1:2 & !is.na(max_stage) & max_stage >= 2),
    peak13  = safe_max(max_stage[day %in% 1:2]),   # 第1~2天最高分期
    .groups = "drop"
  )

grp_levels <- c("Early death (<=72h)", "Early ICU exit", "Transient",
                "Persistent stage 2", "Worsening to stage 3", "Severe from start")

full <- dat %>%
  select(stay_id, heme, tumor, admission_age, sofa_day1, hospital_expire_flag,
         any_rrt, death_28d, death_90d, death_h, t0) %>%
  left_join(d72, by = "stay_id") %>%
  left_join(baseline %>% select(stay_id, gender, charlson_comorbidity_index, sepsis3,
                                metastatic_solid_tumor, renal_disease, congestive_heart_failure,
                                diabetes_without_cc, diabetes_with_cc, creat_baseline_7d),
            by = "stay_id") %>%
  left_join(outc %>% select(stay_id, dod, los_icu, los_hospital), by = "stay_id") %>%
  mutate(
    obs     = ifelse(is.na(obs), 0, obs),
    persist = ifelse(is.na(persist), 0, persist),
    male    = gender %in% c("M", "m", "Male", "male"),
    sepsis3 = as.logical(sepsis3),
    any_rrt = as.logical(any_rrt),
    d28     = as.logical(death_28d),
    d90     = as.logical(death_90d),
    hosp_death = hospital_expire_flag == 1,
    diabetes   = (diabetes_without_cc == 1) | (diabetes_with_cc == 1),
    # 早死：院内死亡有精确时间用 death_h；院外死亡只有日期，日期差<=2天可确定在72小时内
    dod_days0 = as.numeric(as.Date(to_time(dod)) - as.Date(to_time(t0))),
    early72   = (!is.na(death_h) & death_h <= 72) | (is.na(death_h) & !is.na(dod_days0) & dod_days0 <= 2),
    grp_chr = case_when(
      early72                                                          ~ grp_levels[1],
      obs == 0                                                         ~ grp_levels[2],
      persist == 0                                                     ~ grp_levels[3],
      !is.na(peak13) & peak13 >= 3 & (is.na(day0_s) | day0_s < 3)      ~ grp_levels[5],
      !is.na(day0_s) & day0_s >= 3                                     ~ grp_levels[6],
      TRUE                                                             ~ grp_levels[4]
    ),
    grp = factor(grp_chr, levels = grp_levels)
  )

cat("\n各类人数（只用第0~2天划分，界标72小时）:\n")
print(table(full$grp))

# ---------- 2. Table 1 ----------
fmt_med <- function(x) {
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf("%.1f [%.1f-%.1f]", q[2], q[1], q[3])
}
fmt_n <- function(x) sprintf("%d (%.1f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

make_row <- function(label, x, g, type) {
  cells   <- sapply(levels(g), function(l) if (type == "cont") fmt_med(x[g == l]) else fmt_n(x[g == l]))
  overall <- if (type == "cont") fmt_med(x) else fmt_n(x)
  p <- tryCatch({
    if (type == "cont") kruskal.test(x ~ g)$p.value
    else suppressWarnings(chisq.test(table(g, x))$p.value)
  }, error = function(e) NA_real_)
  c(Variable = label, Overall = overall, cells, P = fmt_p(p))
}

g <- full$grp
rows <- list(
  c(Variable = "N", Overall = as.character(nrow(full)), setNames(as.character(as.integer(table(g))), levels(g)), P = ""),
  make_row("Age, years, median [IQR]",            full$admission_age, g, "cont"),
  make_row("Male, n (%)",                         full$male, g, "bin"),
  make_row("Hematologic malignancy, n (%)",       full$heme == 1, g, "bin"),
  make_row("Metastatic solid tumor, n (%)",       full$metastatic_solid_tumor == 1, g, "bin"),
  make_row("Charlson index, median [IQR]",        full$charlson_comorbidity_index, g, "cont"),
  make_row("Congestive heart failure, n (%)",     full$congestive_heart_failure == 1, g, "bin"),
  make_row("Chronic kidney disease, n (%)",       full$renal_disease == 1, g, "bin"),
  make_row("Diabetes, n (%)",                     full$diabetes, g, "bin"),
  make_row("SOFA day 1, median [IQR]",            full$sofa_day1, g, "cont"),
  make_row("Sepsis-3, n (%)",                     full$sepsis3, g, "bin"),
  make_row("Baseline creatinine (7d min), median [IQR]", full$creat_baseline_7d, g, "cont"),
  make_row("KDIGO stage 3 on day 0, n (%)",       full$day0_s >= 3, g, "bin"),
  make_row("Renal replacement therapy, n (%)",    full$any_rrt, g, "bin"),
  make_row("ICU length of stay, days, median [IQR]",      full$los_icu, g, "cont"),
  make_row("Hospital length of stay, days, median [IQR]", full$los_hospital, g, "cont"),
  make_row("Hospital death, n (%)",               full$hosp_death, g, "bin"),
  make_row("28-day death, n (%)",                 full$d28, g, "bin"),
  make_row("90-day death, n (%)",                 full$d90, g, "bin")
)
tab1 <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
cat("\n===== Table 1 =====\n")
print(tab1, row.names = FALSE, right = FALSE)
write.csv(tab1, "table1_by_trajectory.csv", row.names = FALSE)

# ---------- 3. 图1：每一类轨迹、每天的分期构成 ----------
status_lv <- c("Stage 0-1", "Stage 2", "Stage 3", "No data")
cols <- c("#9ecae1", "#f4a21d", "#c0392b", "#e5e5e5")

grid <- expand.grid(stay_id = full$stay_id, day = 0:6) %>%
  left_join(daily, by = c("stay_id", "day")) %>%
  left_join(full %>% select(stay_id, grp), by = "stay_id") %>%
  mutate(status = case_when(
    is.na(max_stage)  ~ "No data",
    max_stage < 2     ~ "Stage 0-1",
    max_stage == 2    ~ "Stage 2",
    TRUE              ~ "Stage 3"
  ))

png("fig1_daily_stage.png", width = 2800, height = 1500, res = 250)
par(mfrow = c(2, 4), mar = c(4, 4, 3, 1))
for (gname in grp_levels) {
  sub  <- grid[grid$grp == gname, ]
  tabg <- table(factor(sub$status, levels = status_lv), factor(sub$day, levels = 0:6))
  prop <- 100 * prop.table(tabg, 2)
  barplot(prop, col = cols, border = NA, names.arg = 0:6, ylim = c(0, 100),
          xlab = "Days since t0", ylab = "% of patients",
          main = paste0(gname, "\n(n = ", length(unique(sub$stay_id)), ")"), cex.main = 0.9)
}
plot.new()
legend("center", legend = status_lv, fill = cols, bty = "n", cex = 1.2,
       title = "Highest KDIGO stage that day\n(creatinine / urine output)")
invisible(dev.off())

# ---------- 4. 图2：72小时界标后的生存曲线（自t0起90天内） ----------
full <- full %>%
  mutate(
    dod_days = as.numeric(as.Date(to_time(dod)) - as.Date(to_time(t0))),
    # 修正：院内死亡有精确时间(death_h)，优先用它；dod 只有日期，用日期相减会把
    # 72~95小时之间死亡的人误当成“3天内死亡”而踢出界标集（3764 vs 应有的约3825）。
    # 院外死亡只有日期：日期差=3天的死亡在48~96小时之间无法确定，按界标之后(3.5天)处理；<=2天已归入早死
    time_d   = ifelse(!is.na(death_h), death_h / 24, ifelse(!is.na(dod_days) & dod_days == 3, 3.5, dod_days)),
    event    = !is.na(time_d) & time_d <= 90
  )

cat("\n===== 生存分析数据自查 =====\n")
cat("事件数（生存分析口径）:", sum(full$event), "  vs  death_90d 为真的人数:", sum(full$d90, na.rm = TRUE), "\n")

lmk <- full %>%
  filter(grp != grp_levels[1], is.na(time_d) | time_d > 3) %>%
  mutate(
    grp  = droplevels(grp),
    tt   = ifelse(event, time_d, 90) - 3
  ) %>%
  filter(tt > 0)

cat("界标集人数:", nrow(lmk), "（应等于 4213 减去早死人数）\n")

fit <- survfit(Surv(tt, event) ~ grp, data = lmk)
lr  <- survdiff(Surv(tt, event) ~ grp, data = lmk)
p_lr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)

cat("\n自t0起第90天的累计生存率（%）：\n")
sm <- summary(fit, times = 87)
print(data.frame(group = levels(lmk$grp), survival_pct = round(100 * sm$surv, 1)), row.names = FALSE)
cat("log-rank P =", fmt_p(p_lr), "\n")

km_cols <- c("#7f7f7f", "#2c7fb8", "#f4a21d", "#d7301f", "#67000d")[seq_along(levels(lmk$grp))]
png("fig2_km_landmark.png", width = 2200, height = 1600, res = 280)
par(mar = c(4.5, 4.5, 3, 1))
plot(fit, col = km_cols, lwd = 2.5, mark.time = FALSE, conf.int = FALSE, ylim = c(0, 1),
     xlab = "Days since landmark (72 h after t0)", ylab = "Survival probability",
     main = "Survival from the 72-hour landmark, by AKI trajectory")
legend("bottomleft", legend = levels(lmk$grp), col = km_cols, lwd = 2.5, bty = "n", cex = 0.9)
text(x = 87 * 0.98, y = 0.97, labels = paste0("Log-rank P ", ifelse(p_lr < 0.001, "< 0.001", paste0("= ", fmt_p(p_lr)))),
     adj = c(1, 1), cex = 0.9)
invisible(dev.off())

cat("\n已保存到：", getwd(), "\n  table1_by_trajectory.csv\n  fig1_daily_stage.png\n  fig2_km_landmark.png\n")
