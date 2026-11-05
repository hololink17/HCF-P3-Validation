# =============================================================
# 10_analysis_16_hierarchical_regression.R —— 分层回归（主口径 N=462）
# -------------------------------------------------------------
# 对应蓝图 R3：α 在协变量之上解释心理症状的增量
# 模型序列（预先指定，不看数据调顺序）：
#   M0  协变量：age_group + gender + education + bmi
#   M1  M0 + α₂（接纳：理论上的近端因子）
#   M2  M0 + 四 α 全量
# 分别对 PHQ-9 与 GAD-7 建模
# 报告：R²、调整 R²、ΔR²、ΔF、ΔF 的 p；连续预测变量 per-SD 标准化以便可比
# 输出：T13_hierarchical_regression.csv、T13b_hierarchical_coefficients.csv
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

L <- load_cohorts()
p <- L[[1]]        # 当前运行队列（Prolific 模式=Prolific；见数模式=见数）

COVARS <- c("age_group", "gender", "education", "bmi")
LABELS <- c("M0 仅协变量", "M1 协变量 + α₂", "M2 协变量 + 四 α 全量")

d <- p %>% select(phq9_total, gad7_total, all_of(ALPHA_AC), all_of(COVARS))
for (cv in COVARS) if (!is.numeric(d[[cv]])) d[[cv]] <- droplevels(factor(d[[cv]]))
n0 <- nrow(d)
d <- d[complete.cases(d), ]
cat(sprintf("\n=== 分层回归（主口径 %d → 完整案例 %d）===\n", n0, nrow(d)))

# ---- per-SD 标准化（α 与 bmi），使系数可比 ----
d <- d %>% mutate(across(all_of(ALPHA_AC), ~ (.x - mean(.x)) / sd(.x)),
                  bmi = (bmi - mean(bmi)) / sd(bmi))

fit_series <- function(dv) {
  rhs0 <- paste(COVARS, collapse = " + ")
  m0 <- lm(as.formula(paste(dv, "~", rhs0)), data = d)
  m1 <- lm(as.formula(paste(dv, "~", rhs0, "+ alpha2")), data = d)
  m2 <- lm(as.formula(paste(dv, "~", rhs0, "+", paste(ALPHA_AC, collapse = " + "))), data = d)
  a01 <- anova(m0, m1)
  a12 <- anova(m1, m2)
  r2  <- c(summary(m0)$r.squared, summary(m1)$r.squared, summary(m2)$r.squared)
  ar2 <- c(summary(m0)$adj.r.squared, summary(m1)$adj.r.squared, summary(m2)$adj.r.squared)
  out <- data.frame(
    dv = dv, model = LABELS,
    R2 = r2, adj_R2 = ar2,
    delta_R2 = c(NA_real_, r2[2] - r2[1], r2[3] - r2[2]),
    delta_F  = c(NA_real_, a01$F[2], a12$F[2]),
    delta_df = c(NA_real_, a01$Df[2], a12$Df[2]),
    delta_p  = c(NA_real_, a01$`Pr(>F)`[2], a12$`Pr(>F)`[2]),
    df_resid = c(m0$df.residual, m1$df.residual, m2$df.residual),
    n = c(nobs(m0), nobs(m1), nobs(m2)), row.names = NULL)
  list(out = out, m2 = m2)
}

res <- list(); models <- list()
for (dv in c("phq9_total", "gad7_total")) {
  r <- fit_series(dv)
  res[[length(res) + 1]] <- r$out
  models[[dv]] <- r$m2
  cat(sprintf("\n  ── 因变量：%s（n=%d）──\n", dv, nobs(r$m2)))
  print(r$out, row.names = FALSE)
}

out <- do.call(rbind, res)
out$delta_R2 <- round(out$delta_R2, 4)
save_tab(out, "T13_hierarchical_regression")

# ---- M2 系数明细（供论文附表）----
coef_rows <- list()
for (dv in names(models)) {
  s <- summary(models[[dv]])$coefficients
  coef_rows[[length(coef_rows) + 1]] <- data.frame(
    dv = dv, term = rownames(s), beta = s[, 1], se = s[, 2],
    t = s[, 3], p = s[, 4], row.names = NULL)
}
save_tab(do.call(rbind, coef_rows), "T13b_hierarchical_coefficients")

cat("\n  ⚠️ 报告纪律：\n")
cat("     · 模型顺序预先指定（协变量 → α₂ → 四 α），未按数据调整\n")
cat("     · α 与 bmi 已 per-SD 标准化，系数可直接比较\n")
cat("     · ΔF 的 p 出自嵌套模型 F 检验，非事后功效\n")
cat("     · 横断面设计：R² 只表示同时关联强度，不作因果解释\n")
