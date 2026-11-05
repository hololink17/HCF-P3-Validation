# =============================================================
# 10_analysis_06_clinical_utility.R —— R4 临床决策效度：DCA 净获益 + NNT
# 对应已发表论文：Condition 4（"经验可验证"→ 决策层面可验证）
# 两队列平行计算，各自按本队列标定切点（内部标定，正文须声明）
# 2026-09-28 修正：α₃ 原始尺度换算公式符号写错（原式把 logit(p)/β₁ 的符号取反）
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts()

net_benefit <- function(y, p, pt) {
  n <- length(y)
  tp <- sum(p >= pt & y == 1); fp <- sum(p >= pt & y == 0)
  tp/n - fp/n * (pt / (1 - pt))
}
dca_table <- function(y, p, pts = seq(0.05, 0.60, by = 0.05)) {
  data.frame(pt = pts,
             model = sapply(pts, function(t) net_benefit(y, p, t)),
             treat_all = sapply(pts, function(t) mean(y) - (1-mean(y)) * (t/(1-t))),
             treat_none = 0)
}
nnt_at <- function(y, p, pt) {
  sel <- p >= pt
  if (sum(sel) == 0) return(c(nnt = NA, n_sel = 0, arr = NA))
  arr <- mean(y[sel]) - mean(y[!sel])
  c(nnt = 1/abs(arr), n_sel = sum(sel), arr = arr)
}

rows <- list()
for (nm in names(L)) {
  d <- L[[nm]] %>% filter(is.finite(alpha3), is.finite(ibs_positive))
  m <- suppressWarnings(glm(ibs_positive ~ alpha3, data = d, family = binomial))
  p <- predict(m, type = "response")
  b0 <- unname(coef(m)[1]); b1 <- unname(coef(m)[2])
  r <- pROC::roc(d$ibs_positive, p, quiet = TRUE)
  pt_best <- as.numeric(pROC::coords(r, "best", ret = "threshold")[1,1])
  # logit(p) = b0 + b1·α₃  →  α₃ = (logit(p) − b0) / b1
  a3_at_pt <- (log(pt_best / (1 - pt_best)) - b0) / b1
  nb <- dca_table(d$ibs_positive, p)
  cat(sprintf("\n=== %s：DCA / NNT（n=%d，IBS+=%d）===\n", nm, nrow(d), sum(d$ibs_positive)))
  cat(sprintf("  α₃ 逻辑回归: β₀=%.3f  β₁=%.3f（β₁<0 = α₃ 越高、IBS 风险越低）\n", b0, b1))
  cat(sprintf("  Youden 最优切点 p = %.4f  →  对应 α₃ 原始尺度 = %.4f\n", pt_best, a3_at_pt))
  cat(sprintf("  预测概率范围: %.3f – %.3f（低于某阈值则无人被选中）\n", min(p), max(p)))
  print(round(nb, 4))
  for (t in c(0.15, 0.20, 0.25, 0.30, 0.40)) {
    v <- nnt_at(d$ibs_positive, p, t)
    cat(sprintf("  阈值 %.2f → NNT=%-6s 选中 n=%-4.0f ARR=%s\n", t,
                ifelse(is.na(v["nnt"]), "NA", sprintf("%.1f", v["nnt"])), v["n_sel"],
                ifelse(is.na(v["arr"]), "NA", sprintf("%.3f", v["arr"]))))
    rows[[length(rows)+1]] <- data.frame(cohort=nm, threshold=t, NNT=v["nnt"],
                                         n_selected=v["n_sel"], ARR=v["arr"],
                                         youden_pt=pt_best, a3_at_youden=a3_at_pt)
  }
  save_tab(cbind(cohort=nm, nb), paste0("T4_dca_", nm))
}
save_tab(do.call(rbind, rows), "T4_nnt_summary")
