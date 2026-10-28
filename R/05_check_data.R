# ============================================================
# 05  数据体检脚本（在你自己电脑的 R / RStudio 里运行）
# 作用：确认四个 CSV 导出正确。只输出汇总数字，不含任何病人记录，
#       Prints summary numbers only (no patient records).
# 使用：1) 把四个 CSV 放进同一个文件夹
#       2) 把下面 data_dir 改成那个文件夹的路径（Windows 用 / 不用 \）
#       3) 全选运行
# ============================================================

# 第一次使用需要安装（只装一次）：
# install.packages(c("readr", "dplyr"))

library(readr)
library(dplyr)

data_dir <- "~/Documents/aki_cancer"   # <<< 改成你的文件夹

baseline <- read_csv(file.path(data_dir, "baseline.csv"),         show_col_types = FALSE)
creat    <- read_csv(file.path(data_dir, "creatinine_raw.csv"),   show_col_types = FALSE)
kdigo    <- read_csv(file.path(data_dir, "kdigo_series.csv"),     show_col_types = FALSE)
outc     <- read_csv(file.path(data_dir, "outcomes.csv"),         show_col_types = FALSE)

cat("\n===== 1. 行数 / 人数 =====\n")
cat("baseline 行数:", nrow(baseline), " 不重复 stay_id:", n_distinct(baseline$stay_id), "  (应为 4213 / 4213)\n")
cat("outcomes 行数:", nrow(outc),     " 不重复 stay_id:", n_distinct(outc$stay_id),     "  (应为 4213 / 4213)\n")
cat("creatinine_raw 行数:", nrow(creat), " 覆盖人数:", n_distinct(creat$stay_id), "\n")
cat("kdigo_series  行数:", nrow(kdigo), " 覆盖人数:", n_distinct(kdigo$stay_id), "\n")

cat("\n===== 2. 肌酐>=3次的人数（应为 2151） =====\n")
print(table(baseline$ge3_cr, useNA = "ifany"))

cat("\n===== 3. 肿瘤类型 =====\n")
cat("含血液肿瘤:", sum(baseline$heme == 1), "  仅实体瘤:", sum(baseline$solid == 1 & baseline$heme == 0), "\n")

cat("\n===== 4. 基线表缺失情况（缺失比例，%） =====\n")
miss <- sapply(baseline, function(x) round(mean(is.na(x)) * 100, 1))
print(miss[miss > 0])

cat("\n===== 5. 结局 =====\n")
cat("住院死亡率 %:", round(mean(outc$hospital_expire_flag == 1) * 100, 1), "\n")
cat("接受透析 %:",   round(mean(outc$any_rrt == TRUE) * 100, 1), "\n")
cat("28天死亡 %:",   round(mean(outc$death_28d == TRUE) * 100, 1), "\n")
cat("90天死亡 %:",   round(mean(outc$death_90d == TRUE) * 100, 1), "\n")

cat("\n===== 6. 时间窗是否正确 =====\n")
cat("creatinine_raw hours_from_t0 范围:", paste(round(range(creat$hours_from_t0), 1), collapse = " ~ "), " (应在 -48 ~ 168 之内)\n")
cat("kdigo_series  hours_from_t0 范围:", paste(round(range(kdigo$hours_from_t0), 1), collapse = " ~ "), " (应在 -24 ~ 168 之内)\n")

cat("\n===== 7. kdigo_series 实际有哪些列 =====\n")
print(names(kdigo))

cat("\n===== 8. 每人 kdigo 行数的分布（决定能不能做轨迹分型） =====\n")
print(summary(as.numeric(table(kdigo$stay_id))))

cat("\n===== 9. 基线肌酐缺失人数 =====\n")
cat(sum(is.na(baseline$creat_baseline_7d)), "人没有 t0 前 7 天的肌酐\n")

cat("\n===== 10. 基线表的列名 =====\n")
print(names(baseline))
