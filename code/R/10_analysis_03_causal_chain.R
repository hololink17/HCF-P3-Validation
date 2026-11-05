# =============================================================
# 10_analysis_03_causal_chain.R —— α₁→α₂→α₃→α₄ 的下游链式回归
# 仅使用后台直出读数进行通用统计，不重构任何内部计分。
# 两队列均适用
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts()
for (nm in names(L)) {
  d <- L[[nm]]; cat(sprintf("\n=== %s：因果链检验 ===\n", nm))
  # 1) 逐段回归：α_k ~ α_{k-1}
  for (k in 2:4) {
    f <- lm(as.formula(paste0("alpha", k, " ~ alpha", k-1)), data = d)
    b <- summary(f)$coefficients
    cat(sprintf("  α%d → α%d : β=%.3f  SE=%.3f  t=%.1f  p=%.2e  R²=%.3f\n",
                k-1, k, b[2,1], b[2,2], b[2,3], b[2,4], summary(f)$r.squared))
  }
}
