
# =============================================================================
# section15.R
# 第 15 节:癌症内部的峰值分期与轨迹分析(v3 稿)
# 运行顺序:先运行 reproduce_all_recovered.R,再运行本文件。
# 需要上一步留下的对象:grp、kid、kdn、cov_names、chk、chk_or、mk、RES、OUT_DIR、CHK
# 预期:最后打印 "OK 108",没有 DIFF 行。
# =============================================================================
suppressPackageStartupMessages({library(dplyr); library(lmtest); library(sandwich); library(pROC)})
gp <- grp
gp$g <- relevel(factor(gp$g), ref = "A_2_transient")
ca <- subset(gp, cancer == 1)
chk("15 癌症 n", nrow(ca), 2656, tol = 0)
G5 <- levels(gp$g)

# (a) 粗死亡率
s15 <- ca |> group_by(g) |> summarise(n = n(), d90 = 100 * mean(d90), .groups = "drop") |> as.data.frame()
rownames(s15) <- s15$g
for (i in seq_along(G5)) {
  chk(paste("15 n", G5[i]),    s15[G5[i], "n"],   c(891, 98, 866, 235, 566)[i], tol = 0)
  chk(paste("15 d90%", G5[i]), s15[G5[i], "d90"], c(31.0, 41.8, 41.6, 54.9, 54.1)[i], tol = 0.051)
}

# (b) 模型、似然比检验、OR
f15 <- list(m0 = reformulate(cov_names, "d90"),
            m1 = reformulate(c(cov_names, "peak3"), "d90"),
            m2 = reformulate(c(cov_names, "g"), "d90"))
g1 <- glm(f15$m1, ca, family = binomial); g2 <- glm(f15$m2, ca, family = binomial)
lrt <- anova(g1, g2, test = "LRT")
chk("15 LRT chisq", lrt$Deviance[2], 20.658, tol = 0.01)
chk("15 LRT df", lrt$Df[2], 3, tol = 0)
chk("15 LRT P", lrt$`Pr(>Chi)`[2], 0.000124, tol = 0.000006)
o15 <- cbind(exp(coef(g2)), exp(confint.default(g2))); colnames(o15) <- c("OR", "lo", "hi")
rn <- c("gB_3_resolved", "gC_2_persistent", "gD_3_worsening", "gE_3_severe")
chk_or("15 OR vs 峰值2期一过性", o15,
       mk(rn, c(1.47, 1.53, 2.28, 1.99), c(0.95, 1.25, 1.69, 1.58), c(2.27, 1.87, 3.09, 2.50)))

# (c) 交叉验证 AUC(沿用 grp$fold)
cv15 <- function(f, d) { p <- rep(NA_real_, nrow(d))
  for (k in sort(unique(d$fold))) {
    fit <- glm(f, d[d$fold != k, ], family = binomial)
    p[d$fold == k] <- predict(fit, d[d$fold == k, ], type = "response") }
  p }
P15 <- lapply(f15, cv15, d = ca)
au <- sapply(P15, function(p) as.numeric(auc(ca$d90, p, quiet = TRUE)))
chk("15 AUC 仅协变量", au["m0"], 0.625, tol = 0.0006)
chk("15 AUC +峰值分期", au["m1"], 0.638, tol = 0.0006)
chk("15 AUC +峰值分期+轨迹", au["m2"], 0.646, tol = 0.0006)
set.seed(2)
dA <- replicate(2000, { i <- sample(nrow(ca), replace = TRUE)
  as.numeric(auc(ca$d90[i], P15$m2[i], quiet = TRUE)) - as.numeric(auc(ca$d90[i], P15$m1[i], quiet = TRUE)) })
ci15 <- quantile(dA, c(.025, .975))
chk("15 AUC 差值", au["m2"] - au["m1"], 0.0083, tol = 0.0006)
chk("15 AUC 差值 lo", ci15[1], -0.0002, tol = 0.0006)
chk("15 AUC 差值 hi", ci15[2], 0.0162, tol = 0.0006)

# (d) 两个峰值层内:未缓解 vs 缓解
ca$unres <- as.integer(ca$g %in% c("C_2_persistent", "D_3_worsening", "E_3_severe"))
mm <- glm(reformulate(c(cov_names, "peak3", "unres", "peak3:unres"), "d90"), ca, family = binomial)
bb <- coef(mm); VV <- vcov(mm)
lc <- function(w) { est <- sum(bb[names(w)] * w); se <- sqrt(drop(w %*% VV[names(w), names(w)] %*% w))
  c(OR = exp(est), lo = exp(est - 1.96 * se), hi = exp(est + 1.96 * se)) }
l2 <- lc(c(unres = 1)); l3 <- lc(c(unres = 1, "peak3:unres" = 1))
chk_or("15 层内未缓解 峰值2期", rbind(a = l2), mk("a", 1.53, 1.25, 1.87))
chk_or("15 层内未缓解 峰值3期", rbind(a = l3), mk("a", 1.41, 0.91, 2.19))
chk("15 层内 交互 P", summary(mm)$coef["peak3:unres", "Pr(>|z|)"], 0.744, tol = 0.0006)

# (e) 合并样本:癌症 x 分组交互(聚类稳健 Wald)
mA <- glm(reformulate(c(cov_names, "g", "cancer"), "d90"), gp, family = binomial)
mB <- glm(reformulate(c(cov_names, "g", "cancer", "cancer:g"), "d90"), gp, family = binomial)
w15 <- waldtest(mA, mB, vcov = function(x) vcovCL(x, cluster = gp$subclass))
chk("15 癌症×分组 交互 F", w15$F[2], 2.4304, tol = 0.001)
chk("15 癌症×分组 交互 P", w15$`Pr(>F)`[2], 0.0455, tol = 0.0006)

# (f) 癌症内部的肾脏结局
gk15 <- left_join(gp, kid, by = "stay_id")
cak <- subset(gk15, cancer == 1)
cr <- cak |> group_by(g) |>
  summarise(rrt_late = 100 * mean(y_rrt_late, na.rm = TRUE), nonrec = 100 * mean(y_nonrec, na.rm = TRUE),
            make = 100 * mean(y_make), .groups = "drop") |> as.data.frame()
rownames(cr) <- cr$g
ep <- list(rrt_late = c(0.9, 2.1, 1.7, 4.6, 5.3), nonrec = c(6.9, 9.2, 12.0, 23.6, 31.0),
           make = c(15.0, 18.4, 24.8, 44.3, 47.5))
for (k in names(ep)) for (i in seq_along(G5))
  chk(paste("15 癌症内", k, "%", G5[i]), cr[G5[i], k], ep[[k]][i], tol = 0.051)
or15 <- function(y) { d <- cak[!is.na(cak[[y]]), ]
  m <- glm(reformulate(c(cov_names, "g"), y), d, family = binomial)
  tb <- cbind(exp(coef(m)), exp(confint.default(m))); colnames(tb) <- c("OR", "lo", "hi"); tb }
chk_or("15 癌症内 MAKE", or15("y_make"),
       mk(rn, c(1.16, 1.88, 3.67, 3.40), c(0.65, 1.46, 2.62, 2.62), c(2.04, 2.41, 5.13, 4.43)))
chk_or("15 癌症内 RRT>72h", or15("y_rrt_late"),
       mk(rn, c(2.03, 2.03, 3.91, 3.78), c(0.41, 0.85, 1.47, 1.65), c(9.98, 4.86, 10.38, 8.62)))
chk_or("15 癌症内 未恢复", or15("y_nonrec"),
       mk(rn, c(1.42, 1.89, 3.84, 4.45), c(0.61, 1.28, 2.31, 3.02), c(3.30, 2.78, 6.38, 6.56)))
chk_or("15 癌症内 任何RRT", or15("y_rrt_any"),
       mk(rn, c(5.40, 2.40, 9.64, 9.87), c(1.81, 1.08, 4.36, 4.82), c(16.09, 5.35, 21.32, 20.21)))

# (g) 排除 hospice 后未恢复的 n/N(需要第 12 节生成的 kdn)
hn <- kdn |> group_by(cancer) |>
  summarise(n = sum(!is.na(y_nonrec)), ev = sum(y_nonrec, na.rm = TRUE), .groups = "drop")
chk("15 排除hospice 未恢复 n 肿瘤", hn$n[hn$cancer == 1], 1722, tol = 0)
chk("15 排除hospice 未恢复 事件 肿瘤", hn$ev[hn$cancer == 1], 227, tol = 0)
chk("15 排除hospice 未恢复 n 非肿瘤", hn$n[hn$cancer == 0], 4062, tol = 0)
chk("15 排除hospice 未恢复 事件 非肿瘤", hn$ev[hn$cancer == 0], 497, tol = 0)

# 输出
RES$sec15_crude <- s15
RES$sec15_or    <- as.data.frame(round(o15, 3))
RES$sec15_auc   <- data.frame(model = names(au), auc = round(au, 4))
write.csv(CHK, file.path(OUT_DIR, "checks_with_sec15.csv"), row.names = FALSE)
s15chk <- CHK[grepl("^15 ", CHK$item), ]
print(table(s15chk$status))
print(s15chk[s15chk$status == "DIFF", ], row.names = FALSE)

