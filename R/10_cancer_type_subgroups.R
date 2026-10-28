# ============================================================
# 10  按肿瘤类型的亚组分析（新增）
# 运行顺序：先跑完 07（会留下 md_unified、md_surv、cl_or、cov_str、lev72、p72），再跑本脚本。
#
# 亚组定义（互斥，按优先级）：
#   Hematologic            : dat_c$heme == 1
#   Metastatic solid tumor : 非血液肿瘤，且 Charlson 的 metastatic_solid_tumor == 1
#   Non-metastatic solid   : 其余癌症患者
# 匹配后每个癌症患者有 2 个对照（同一 subclass），对照继承其所属癌症患者的亚组标签，
# 所以每个亚组内部仍是“匹配后的癌症 vs 对照”比较。
#
# 分析（与主分析一致：72小时界标，早死已剔除，聚类稳健SE按 subclass）：
#   A 亚组构成和各亚组的轨迹分布
#   B 各亚组内 癌症 vs 非癌症 的 OR（院内死亡、90天死亡）
#   C 亚组异质性：cancer × 亚组 交互项 Wald P
#   D 各亚组 × 合并轨迹（早期退出 / 短暂 / 持续2期 / 3期型）的 90 天 OR
#   E 各亚组粗死亡率
# Prints summary numbers only; also saves 4 summary csv files.
# ============================================================

library(dplyr)
library(sandwich)
library(lmtest)

need <- c("md_unified", "md_surv", "cl_or", "cov_str", "dat_c")
miss <- need[!sapply(need, exists)]
if (length(miss) > 0) stop("缺少对象：", paste(miss, collapse = ", "), "。请先完整运行 07 脚本。")

# ---------- 0. 看一眼 dat_c 里的肿瘤变量（只是汇总计数） ----------
cat("\n===== 0. 肿瘤相关变量检查 =====\n")
cat("dat_c 列名中含 heme/tumor/cancer 的：", paste(grep("heme|tumor|cancer|onc", names(dat_c), value = TRUE), collapse = ", "), "\n")
if ("tumor" %in% names(dat_c)) { cat("table(dat_c$tumor):\n"); print(head(sort(table(dat_c$tumor, useNA = "ifany"), decreasing = TRUE), 15)) }
cat("table(dat_c$heme):\n"); print(table(dat_c$heme, useNA = "ifany"))

# 转移性实体瘤标记：优先用 baseline_c，其次 baseline，再看 dat_c 里有没有
bl <- NULL
for (nm in c("baseline_c", "baseline")) if (exists(nm)) { bl <- get(nm); break }
met <- NULL
if (!is.null(bl) && "metastatic_solid_tumor" %in% names(bl)) {
  met <- bl %>% select(stay_id, metastatic_solid_tumor)
} else if ("metastatic_solid_tumor" %in% names(dat_c)) {
  met <- dat_c %>% select(stay_id, metastatic_solid_tumor)
} else stop("Column metastatic_solid_tumor not found in baseline_c, baseline or dat_c.")

# ---------- 1. 构造亚组标签（只对癌症患者），再传给匹配对照 ----------
sub_c <- dat_c %>% select(stay_id, heme) %>%
  left_join(met, by = "stay_id") %>%
  mutate(subtype = case_when(
    heme == 1                       ~ "Hematologic",
    metastatic_solid_tumor == 1     ~ "Metastatic solid",
    TRUE                            ~ "Non-metastatic solid")) %>%
  select(stay_id, subtype)

lab_ids <- md_unified %>% filter(cancer == 1) %>% select(stay_id, subclass) %>%
  left_join(sub_c, by = "stay_id") %>% select(subclass, subtype) %>% distinct()
if (any(duplicated(lab_ids$subclass))) stop("同一个 subclass 里出现多个癌症患者，说明不是 1:2 匹配，请告诉我。")
if (any(is.na(lab_ids$subtype)))        stop("有癌症患者没有亚组标签，请检查 dat_c 是否包含全部癌症 stay_id。")

sub_levels <- c("Hematologic", "Non-metastatic solid", "Metastatic solid")
add_sub <- function(d) {
  d %>% select(-any_of("subtype")) %>% left_join(lab_ids, by = "subclass") %>%
    mutate(subtype = factor(subtype, levels = sub_levels))
}
mu <- add_sub(md_unified)
ms <- add_sub(md_surv)
ms$pheno72 <- droplevels(ms$pheno72)

# 合并轨迹（亚组内样本量有限，把恶化型和起病即重度合并）
comb <- function(x) factor(case_when(
  x == "2_no_data_day1to2"                                  ~ "Early ICU exit",
  x == "3_transient"                                        ~ "Transient",
  x == "5_persistent_stage2"                                ~ "Persistent stage 2",
  x %in% c("6_worsening_to_stage3", "7_severe_from_start")  ~ "Stage 3 (worsening or from start)"),
  levels = c("Early ICU exit", "Transient", "Persistent stage 2", "Stage 3 (worsening or from start)"))
ms$traj4 <- comb(as.character(ms$pheno72))

# ---------- A. 构成与轨迹分布 ----------
cat("\n===== A. 癌症患者亚组构成（匹配队列，含早死） =====\n")
tc <- table(mu$subtype[mu$cancer == 1])
print(tc); print(round(100 * prop.table(tc), 1))
cat("\n癌症患者：亚组 × 轨迹（行 %，匹配队列含早死）\n")
tt <- table(mu$subtype[mu$cancer == 1], mu$pheno72[mu$cancer == 1])
print(tt); print(round(100 * prop.table(tt, 1), 1))
chA <- suppressWarnings(chisq.test(tt)); cat("亚组间轨迹分布差异：卡方 =", round(chA$statistic, 1), " df =", chA$parameter, " P =", signif(chA$p.value, 3), "\n")
write.csv(as.data.frame.matrix(tt), "subgroup_traj_counts.csv")

# ---------- B. 各亚组内 癌症 vs 非癌症 ----------
cat("\n===== B. 各亚组：癌症 vs 非癌症 的调整 OR（界标后，早死已剔除） =====\n")
rowB <- lapply(sub_levels, function(s) {
  d <- ms[ms$subtype == s, ]
  out <- data.frame(Subgroup = s, N_Cancer = sum(d$cancer == 1), N_Control = sum(d$cancer == 0))
  for (y in c("hosp_death", "d90")) {
    f <- tryCatch(glm(as.formula(paste(y, "~ cancer +", cov_str)), data = d, family = binomial), error = function(e) NULL)
    r <- if (is.null(f)) c(NA, NA, NA) else cl_or(f, d)
    out[[paste0(ifelse(y == "d90", "D90", "Hosp"), c("_OR", "_Lo", "_Hi"))[1]]] <- r[1]
    out[[paste0(ifelse(y == "d90", "D90", "Hosp"), "_Lo")]] <- r[2]
    out[[paste0(ifelse(y == "d90", "D90", "Hosp"), "_Hi")]] <- r[3]
  }
  out
})
resB <- do.call(rbind, rowB); print(resB)
write.csv(resB, "subgroup_or_by_cancer_type.csv", row.names = FALSE)

# ---------- C. 异质性：cancer × 亚组 ----------
cat("\n===== C. 异质性检验：cancer × 亚组 交互项（Wald, 聚类稳健） =====\n")
for (y in c("hosp_death", "d90")) {
  fa <- glm(as.formula(paste(y, "~ cancer + subtype +", cov_str)), data = ms, family = binomial)
  fi <- glm(as.formula(paste(y, "~ cancer * subtype +", cov_str)), data = ms, family = binomial)
  pv <- waldtest(fa, fi, test = "Chisq", vcov = function(x) vcovCL(x, cluster = ms$subclass))$`Pr(>Chisq)`[2]
  cat(y, ": P for heterogeneity =", signif(pv, 3), "\n")
}

# ---------- D. 亚组 × 合并轨迹 ----------
cat("\n===== D. 各亚组 × 合并轨迹：癌症 vs 非癌症 的 90 天 OR =====\n")
rowD <- list()
for (s in sub_levels) for (tr in levels(ms$traj4)) {
  d <- ms[ms$subtype == s & ms$traj4 == tr, ]
  nc <- sum(d$cancer == 1); nk <- sum(d$cancer == 0)
  r <- c(NA, NA, NA)
  if (nc >= 20 && nk >= 20 && sum(d$d90) >= 10) {
    f <- tryCatch(suppressWarnings(glm(as.formula(paste("d90 ~ cancer +", cov_str)), data = d, family = binomial)), error = function(e) NULL)
    if (!is.null(f) && f$converged) r <- cl_or(f, d)
  }
  rowD[[length(rowD) + 1]] <- data.frame(Subgroup = s, Trajectory = tr, N_Cancer = nc, N_Control = nk,
                                         D90_OR = r[1], D90_Lo = r[2], D90_Hi = r[3])
}
resD <- do.call(rbind, rowD); print(resD)
cat("（样本太少或不收敛的格子显示 NA；n<20 或死亡<10 不拟合）\n")
write.csv(resD, "subgroup_by_trajectory_or.csv", row.names = FALSE)

# ---------- E. 粗死亡率 ----------
cat("\n===== E. 各亚组粗死亡率 (%)（界标后） =====\n")
crudeE <- ms %>% group_by(subtype, cancer) %>%
  summarise(n = n(), hosp = round(100 * mean(hosp_death), 1), d90 = round(100 * mean(d90), 1), .groups = "drop") %>%
  as.data.frame()
print(crudeE)
write.csv(crudeE, "subgroup_crude_mortality.csv", row.names = FALSE)

cat("\n===== Done =====\n")
