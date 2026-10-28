# ============================================================
# 06  规则法 AKI 轨迹分型（第一版）
# 前提：R 里已经有 baseline、outc、kdigo 三个表（上一步读好的），
#       并且已经运行过 as.logical 那四行。
# Output: prints class counts and outcome proportions (summary numbers only).
# Patient-level results are saved in the working directory only; do not upload them.
# ============================================================

library(dplyr)

# ---------- 0. 小工具 ----------
safe_max <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA_real_ else max(x) }
to_time  <- function(x) if (inherits(x, "POSIXt")) x else as.POSIXct(as.character(x), tz = "UTC")

# ---------- 1. 每小时 KDIGO 分期 -> 每天最高分期（第0~6天） ----------
# 说明：分期只用肌酐和尿量，不含透析（和定义 t0 时一致）。
#       NA 表示该行没有测量，0 表示测了但没有 AKI。
daily <- kdigo %>%
  filter(hours_from_t0 >= 0) %>%
  mutate(
    stage_h = pmax(as.numeric(aki_stage_creat), as.numeric(aki_stage_uo), na.rm = TRUE),
    day     = pmin(floor(hours_from_t0 / 24), 6)          # 168小时并入第6天
  ) %>%
  group_by(stay_id, day) %>%
  summarise(max_stage = safe_max(stage_h), .groups = "drop")

# ---------- 2. 每人一行的摘要特征 ----------
traj <- daily %>%
  group_by(stay_id) %>%
  summarise(
    day0_stage  = safe_max(max_stage[day == 0]),
    obs_1_3     = sum(day %in% 1:3 & !is.na(max_stage)),                 # 第1~3天有多少天有数据
    persist_1_3 = sum(day %in% 1:3 & !is.na(max_stage) & max_stage >= 2),# 其中有多少天仍>=2期
    peak_later  = safe_max(max_stage[day >= 1]),                         # 第1~6天最高分期
    late_aki    = any(day %in% 4:6 & !is.na(max_stage) & max_stage >= 2),# 第4~6天是否再次>=2期
    .groups = "drop"
  )

# ---------- 3. 合并基线和结局 ----------
dat <- baseline %>%
  select(stay_id, heme, solid, admission_age, sofa_day1) %>%
  left_join(
    outc %>% select(stay_id, t0, deathtime, hospital_expire_flag, any_rrt, death_28d, death_90d),
    by = "stay_id"
  ) %>%
  left_join(traj, by = "stay_id") %>%
  mutate(
    obs_1_3     = ifelse(is.na(obs_1_3), 0, obs_1_3),
    persist_1_3 = ifelse(is.na(persist_1_3), 0, persist_1_3),
    late_aki    = ifelse(is.na(late_aki), FALSE, late_aki),
    death_h     = as.numeric(difftime(to_time(deathtime), to_time(t0), units = "hours"))
  )

# ---------- 4. 规则分型（按顺序判断，先满足的先归类） ----------
dat <- dat %>%
  mutate(
    pheno = case_when(
      !is.na(death_h) & death_h <= 72                                   ~ "1_early_death_72h",
      obs_1_3 == 0                                                      ~ "2_no_data_day1to3",
      persist_1_3 == 0 & !late_aki                                      ~ "3_transient",
      persist_1_3 == 0 &  late_aki                                      ~ "4_recurrent_late",
      !is.na(peak_later) & peak_later >= 3 & (is.na(day0_stage) | day0_stage < 3) ~ "6_worsening_to_stage3",
      !is.na(day0_stage) & day0_stage >= 3                              ~ "7_severe_from_start",
      TRUE                                                              ~ "5_persistent_stage2"
    ),
    tumor = ifelse(heme == 1, "heme", "solid_only")
  )

# ---------- 5. 输出汇总 ----------
cat("\n===== A. 第0天分期分布（t0当天最高分期） =====\n")
print(table(dat$day0_stage, useNA = "ifany"))

cat("\n===== B. 各类轨迹的人数和结局 =====\n")
res <- dat %>%
  group_by(pheno) %>%
  summarise(
    n          = n(),
    pct        = round(100 * n() / nrow(dat), 1),
    hosp_death = round(100 * mean(hospital_expire_flag == 1, na.rm = TRUE), 1),
    death_28d  = round(100 * mean(death_28d, na.rm = TRUE), 1),
    death_90d  = round(100 * mean(death_90d, na.rm = TRUE), 1),
    rrt        = round(100 * mean(any_rrt, na.rm = TRUE), 1),
    heme_pct   = round(100 * mean(heme == 1, na.rm = TRUE), 1),
    sofa_med   = median(sofa_day1, na.rm = TRUE),
    age_med    = median(admission_age, na.rm = TRUE),
    .groups = "drop"
  )
print(as.data.frame(res))
cat("各类人数合计:", sum(res$n), " (应为 4213)\n")

cat("\n===== C. 血液肿瘤 vs 仅实体瘤：各类轨迹所占比例（%） =====\n")
print(round(100 * prop.table(table(dat$pheno, dat$tumor), margin = 2), 1))

cat("\n===== D. 缺失检查 =====\n")
cat("死亡时间 deathtime 非空的人数:", sum(!is.na(dat$deathtime)),
    " 住院死亡人数:", sum(dat$hospital_expire_flag == 1, na.rm = TRUE), "\n")
cat("没有任何 >=0 小时 KDIGO 行的人数:", sum(!dat$stay_id %in% daily$stay_id), "\n")

# ---------- 6. Save locally (do not upload) ----------
saveRDS(daily, "daily_stage.rds")
saveRDS(dat,   "traj_rule_based.rds")
cat("\n已保存 daily_stage.rds 和 traj_rule_based.rds 到：", getwd(), "\n")
