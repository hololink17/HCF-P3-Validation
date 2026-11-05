# =============================================================
# 10_analysis_02_alpha_structure.R —— R2 构念结构：α 相关矩阵 + VIF + 锚定偏差
# 对应已发表论文：Condition 1/2（公理与状态定义）；回应"α 只是总分冗余"的质疑
# 两队列均适用
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts(); AC <- c("alpha1","alpha2","alpha3","alpha4")
for (nm in names(L)) {
  d <- L[[nm]]
  cat(sprintf("\n=== %s：α 相关矩阵 ===\n", nm))
  cm <- cor(d[, AC], use = "pairwise.complete.obs"); print(round(cm, 3))
  cat("  VIF: ")
  fit <- lm(alpha3 ~ alpha1 + alpha2 + alpha4, data = d)
  v <- 1/(1 - summary(fit)$r.squared)
  cat(sprintf("alpha3 对其它三α 的 VIF = %.2f\n", v))
  # 锚定偏差：算法 α vs 主观自评
  if (all(c("anchor1","anchor2","anchor3","anchor4") %in% names(d))) {
    cat("  锚定偏差（算法α vs 主观自评 r）: ")
    for (i in 1:4) {
      a <- d[[paste0("alpha", i)]]; b <- d[[paste0("anchor", i)]]
      cat(sprintf("α%d=%.3f  ", i, cor(a, b, use = "pairwise.complete.obs")))
    }
    cat("\n")
  }
  # 单因子解释力：第一主成分方差占比（检验"一个总因子是否解释一切"）
  pc <- prcomp(na.omit(d[, AC]), scale. = TRUE)
  cat(sprintf("  PCA 第一主成分方差占比 = %.1f%%\n", 100 * pc$sdev[1]^2 / sum(pc$sdev^2)))
}
# 落盘：两队列相关矩阵
out <- do.call(rbind, lapply(names(L), function(nm) {
  d <- L[[nm]]; cm <- cor(d[, AC], use = "pairwise.complete.obs")
  data.frame(cohort = nm, from = rep(rownames(cm), 4), to = rep(colnames(cm), each = 4),
             r = as.vector(cm))
}))
save_tab(out, "T2_alpha_correlation")
