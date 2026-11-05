# =============================================================
# 10_analysis_09_carrier.R —— R7 载体层：f(a)–f(d) 与结局
# 对应已发表论文：Condition 1（公理：双层层级架构 意识 / 载体）
# 两队列均适用（公开基底统一使用 carrier_a–carrier_d）
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts(); CR <- c("carrier_a","carrier_b","carrier_c","carrier_d")
# 注：载体一律以 carrier_a..d 代号呈现，不输出任何中文层级名（对外边界纪律）
for (nm in names(L)) {
  d <- L[[nm]]; cat(sprintf("\n=== %s：载体层（n=%d）===\n", nm, nrow(d)))
  cm <- cor(d[, CR], use = "pairwise.complete.obs"); print(round(cm, 3))
  for (i in seq_along(CR)) {
    r1 <- cor(d[[CR[i]]], d$alpha3, use="complete.obs")
    v  <- auc_oriented(d$ibs_positive, d[[CR[i]]])
    dp <- cohens_d(d[[CR[i]]], d$ibs_positive)
    cat(sprintf("  %-14s  与α₃ r=%+.3f   IBS AUC=%.3f   Cohen's d=%.2f\n",
                CR[i], r1, v["auc"], dp))
  }
  # 载体 vs 意识层：谁的预测力更强
  cv_c <- cv_auc(d, CR, ycol="ibs_positive")
  cv_a <- cv_auc(d, c("alpha1","alpha2","alpha3","alpha4"), ycol="ibs_positive")
  cat(sprintf("  IBS 预测 5折CV：意识层(4α) AUC=%.3f  vs  载体层(4f) AUC=%.3f\n", cv_a["mean"], cv_c["mean"]))
}
