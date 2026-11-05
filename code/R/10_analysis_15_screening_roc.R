# =============================================================
# 10_analysis_15_screening_roc.R —— α 因子筛查抑郁/焦虑的判别力（主口径 N=462）
# -------------------------------------------------------------
# 对应蓝图 R3：筛查判别力（原 9.18 版为 N=466 口径，本模块按 462 重算）
# 设计要点（面向审稿人）：
#   ① 结局切点先验指定：PHQ-9 ≥ 10、GAD-7 ≥ 10（临床常用切点），不看数据选
#   ② 方向先验指定：α 越高 = 状态越健康 → 预测「非病例」，AUC 用 −α 定向，
#      不采用「取 max(auc, 1−auc)」的自动定向（那会掩盖方向错误）
#   ③ 双口径并行报告：表观 AUC（全样本）与 5 折交叉验证 AUC（防乐观偏倚）
#   ④ 不报事后功效（Hoenig & Heisey 2001）
# 输出：T12_screening_roc.csv
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

PHQ_CUT <- 10   # 与 00_shared_rules.R 的 self_ok 同一切点
GAD_CUT <- 10

L <- load_cohorts()
p <- L[[1]]        # 当前运行队列（Prolific 模式=Prolific；见数模式=见数）
cat(sprintf("\n=== α 筛查判别力（主口径 N=%d）===\n", nrow(p)))

p$dep <- as.numeric(p$phq9_total >= PHQ_CUT)
p$anx <- as.numeric(p$gad7_total >= GAD_CUT)
cat(sprintf("  抑郁阳性 PHQ-9≥%d：%d / %d（%.1f%%）\n", PHQ_CUT, sum(p$dep, na.rm = TRUE), nrow(p),
            100 * mean(p$dep, na.rm = TRUE)))
cat(sprintf("  焦虑阳性 GAD-7≥%d：%d / %d（%.1f%%）\n", GAD_CUT, sum(p$anx, na.rm = TRUE), nrow(p),
            100 * mean(p$anx, na.rm = TRUE)))

rows <- list()

# ---------- 单变量 α → 结局 ----------
for (oc in c("dep", "anx")) {
  y <- p[[oc]]
  for (a in ALPHA_AC) {
    r <- auc_ci(y, p[[a]], orient = TRUE)   # 方向先验：−α
    rows[[length(rows) + 1]] <- data.frame(
      outcome = ifelse(oc == "dep", "PHQ-9≥10", "GAD-7≥10"),
      model = a, n_pred = 1,
      auc = r["auc"], lo = r["lo"], hi = r["hi"],
      n = r["n"], pos = r["pos"], method = "表观 AUC（DeLong CI）",
      row.names = NULL)
  }
}

# ---------- 多变量：四 α 联合（logistic）----------
for (oc in c("dep", "anx")) {
  p$.y <- p[[oc]]
  # 表观 AUC
  m <- suppressWarnings(glm(.y ~ alpha1 + alpha2 + alpha3 + alpha4, data = p, family = binomial))
  pr <- suppressWarnings(predict(m, type = "response"))
  ok <- is.finite(pr) & is.finite(p$.y)
  rr <- pROC::roc(p$.y[ok], pr[ok], quiet = TRUE)
  ci <- as.numeric(pROC::ci.auc(rr))
  rows[[length(rows) + 1]] <- data.frame(
    outcome = ifelse(oc == "dep", "PHQ-9≥10", "GAD-7≥10"),
    model = "四 α 联合（20 题）", n_pred = 4,
    auc = as.numeric(pROC::auc(rr)), lo = ci[1], hi = ci[3],
    n = sum(ok), pos = sum(p$.y[ok]), method = "表观 AUC（DeLong CI）",
    row.names = NULL)
  # 交叉验证 AUC（防乐观偏倚）
  cv <- cv_auc(transform(p, ybin = p[[oc]]), preds = ALPHA_AC, ycol = "ybin", k = 5)
  rows[[length(rows) + 1]] <- data.frame(
    outcome = ifelse(oc == "dep", "PHQ-9≥10", "GAD-7≥10"),
    model = "四 α 联合（20 题）", n_pred = 4,
    auc = cv["mean"], lo = cv["mean"] - cv["sd"], hi = cv["mean"] + cv["sd"],
    n = cv["n"], pos = NA_integer_, method = "5 折交叉验证 AUC（mean ± SD）",
    row.names = NULL)
}

out <- do.call(rbind, rows)
save_tab(out, "T12_screening_roc")

# ---------- 落盘日志 ----------
cat("\n  ── 判别力汇总 ──\n")
print(out, row.names = FALSE)
cat("\n  ⚠️ 报告纪律：\n")
cat("     · 表观 AUC 为全样本拟合值，存在乐观偏倚；5 折 CV 是内部验证，可减轻但不消除该偏倚，且不是外部验证\n")
cat("     · 方向为先验指定（−α），未做自动定向\n")
cat("     · 切点 10 为先验指定，非数据驱动\n")
cat("     · 不以同源内部分数预测由相同读数定义的结局，避免循环论证。\n")
