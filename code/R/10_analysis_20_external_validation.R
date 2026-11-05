# =============================================================
# 10_analysis_20_external_validation.R —— R3/R4.5 【零成本补强】跨队列切点外部验证
# -------------------------------------------------------------
# 本模块在独立样本上评价后台输出的下游统计表现：
#   · 除「内部 5 折交叉验证」（模块 15 / 17）外，增加**跨队列外部验证**：
#       训练队列（A）拟合 + 定切点 → 验证队列（B）独立评估 Sens/Spec/PPV/NPV/AUC + 校准
#       再反向做一次（B → A），双向互证
# 为什么这算"外部验证"：Prolific 与见数是
#   不同平台（Prolific 众包 / 见数数据集市）
#   不同采集时间（2026-05～06 / 2026-08～09）
#   不同来源人群 —— 满足外部验证对"独立样本"的基本要求
# 实现口径（2026-09-30 统一）：所有模型（单变量 / 多变量）一律**先概率化**——
#   训练队列拟合 logistic → 得到概率 → Youden 定阈值 → 验证队列外推概率。
#   验证 AUC / CI / 校准全部基于**同一连续概率**，口径自洽。
# 本模块为跨队列分析：无论 P3_COHORT 为哪一种，两侧都写出同一份结果（与 T5/T10/T11 同规）。
# 口径：Prolific N=462（每人一行）｜见数 N=504
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

# ---------- 工具函数 ----------
# 训练队列上按 Youden 指数取概率阈值
youden_cut <- function(y, p) {
  ok <- is.finite(y) & is.finite(p); y <- y[ok]; p <- p[ok]
  if (length(unique(y)) < 2 || length(p) < 20) return(NA_real_)
  r <- tryCatch(pROC::roc(y, p, quiet = TRUE, direction = "<"), error = function(e) NULL)
  if (is.null(r)) return(NA_real_)
  co <- tryCatch(pROC::coords(r, "best", best.method = "youden",
                              ret = c("threshold"), transpose = FALSE),
                 error = function(e) NULL)
  if (is.null(co) || !nrow(co)) return(NA_real_)
  as.numeric(co$threshold[1])
}

# 诊断指标（给定二分类与真值）
diag_metrics <- function(y, pred_bin) {
  ok <- is.finite(y) & !is.na(pred_bin); y <- y[ok]; pred_bin <- pred_bin[ok]
  tp <- sum(pred_bin == 1 & y == 1); fp <- sum(pred_bin == 1 & y == 0)
  fn <- sum(pred_bin == 0 & y == 1); tn <- sum(pred_bin == 0 & y == 0)
  sens <- ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_)
  spec <- ifelse((fp + tn) > 0, tn / (fp + tn), NA_real_)
  ppv  <- ifelse((tp + fp) > 0, tp / (tp + fp), NA_real_)
  npv  <- ifelse((tn + fn) > 0, tn / (tn + fn), NA_real_)
  lrp  <- ifelse(is.na(sens) | is.na(spec) | spec == 1, NA_real_, sens / (1 - spec))
  lrn  <- ifelse(is.na(sens) | is.na(spec) | spec == 0, NA_real_, (1 - sens) / spec)
  c(tp = tp, fp = fp, fn = fn, tn = tn, sens = sens, spec = spec,
    ppv = ppv, npv = npv, lr_pos = lrp, lr_neg = lrn,
    acc = (tp + tn) / length(y))
}

# 校准：calibration slope / intercept（理想 1 / 0）+ Hosmer–Lemeshow（5 组）
# ⚠ 大样本下 HL 检验对微小偏离也高度敏感（n≈500 时几乎必然 p<0.05），
#   判读应以**校准斜率与截距**为主，HL p 仅作参考。
# 🔴 落盘精度（2026-10-06 第二十一轮修）：hl_p 原用 round(…, 4)，会把 <5e-5 的 p 压成 0
#   ——落盘成「0」等于宣称不可能事件，与本文落盘不变量⑤（p 一律先过 %.4g）不符。
#   改用 fmt_p()：≥1e-4 保持 4 位小数，更小者保留 4 位有效数字（科学计数）。
fmt_p <- function(x) {
  if (length(x) != 1 || is.na(x)) return(NA_real_)
  as.numeric(sprintf("%.4g", x))
}
calibration <- function(y, p) {
  ok <- is.finite(y) & is.finite(p) & p > 0 & p < 1
  y <- y[ok]; p <- p[ok]
  if (length(y) < 30 || length(unique(y)) < 2) {
    return(c(cal_slope = NA_real_, cal_intercept = NA_real_, hl_p = NA_real_))
  }
  lg <- log(p / (1 - p))
  m <- tryCatch(suppressWarnings(glm(y ~ lg, family = binomial)), error = function(e) NULL)
  slope <- if (is.null(m)) NA_real_ else unname(coef(m)[2])
  intc  <- if (is.null(m)) NA_real_ else unname(coef(m)[1])
  br <- unique(quantile(p, seq(0, 1, length.out = 6)))
  hl_p <- NA_real_
  if (length(br) == 6) {
    g <- cut(p, br, include.lowest = TRUE)
    obs <- tapply(y, g, sum); exp_ <- tapply(p, g, sum); n_g <- tapply(y, g, length)
    keep <- !is.na(obs)
    if (sum(keep) >= 3) {
      stat <- sum((obs[keep] - exp_[keep])^2 / pmax(exp_[keep], 1e-6)) +
        sum(((n_g[keep] - obs[keep]) - (n_g[keep] - exp_[keep]))^2 / pmax(n_g[keep] - exp_[keep], 1e-6))
      # 🔴 必须用 lower.tail=FALSE（上尾直算）：`1 - pchisq(stat, df)` 在 p 极小时
      #   发生灾难性抵消（低于 ~1e-16 一律归零），会把真实存在的小 p 报成 0。
      hl_p <- pchisq(stat, df = sum(keep) - 2, lower.tail = FALSE)
    }
  }
  c(cal_slope = slope, cal_intercept = intc, hl_p = hl_p)
}

# ---------- 数据准备 ----------
L <- load_both()
for (nm in names(L)) {
  d <- L[[nm]]
  d$phq9_pos <- as.numeric(d$phq9_total >= 10)
  d$gad7_pos <- as.numeric(d$gad7_total >= 10)
  L[[nm]] <- d
}

cat("\n=== R3/R4.5 跨队列切点外部验证 ===\n")
cat(sprintf("  Prolific n=%d（IBS 阳性 %d，%.1f%%）｜ 见数 n=%d（IBS 阳性 %d，%.1f%%）\n",
            nrow(L$prolific), sum(L$prolific$ibs_positive == 1, na.rm = TRUE),
            mean(L$prolific$ibs_positive == 1, na.rm = TRUE) * 100,
            nrow(L$jianshu),  sum(L$jianshu$ibs_positive == 1, na.rm = TRUE),
            mean(L$jianshu$ibs_positive == 1, na.rm = TRUE) * 100))
cat("  方法：训练队列拟合 logistic → Youden 定概率阈值 → 验证队列外推（AUC/CI/校准同源）\n")

OUTCOME_SPECS <- list(
  list(outcome = "IBS 阳性", ycol = "ibs_positive", preds = "alpha3"),
  list(outcome = "PHQ-9≥10", ycol = "phq9_pos",     preds = c("alpha1", "alpha2", "alpha3", "alpha4")),
  list(outcome = "GAD-7≥10", ycol = "gad7_pos",     preds = c("alpha1", "alpha2", "alpha3", "alpha4"))
)
DIRS <- list(c("prolific", "jianshu"), c("jianshu", "prolific"))

rows <- list()
for (dr in DIRS) {
  tr_name <- dr[1]; te_name <- dr[2]
  tr <- L[[tr_name]]; te <- L[[te_name]]
  for (o in OUTCOME_SPECS) {
    models <- list(list(preds = o$preds[1], label = sprintf("单变量 %s", o$preds[1])))
    if (length(o$preds) > 1)
      models <- c(models, list(list(preds = o$preds, label = "四 α 联合")))

    for (mm in models) {
      preds <- mm$preds
      # 方向对齐：使预测子数值越高 = 风险越高（方向由**训练队列**判定，避免信息泄漏）
      cor_tr <- suppressWarnings(cor(tr[[preds[1]]], tr[[o$ycol]], use = "complete.obs"))
      sgn <- if (!is.na(cor_tr) && cor_tr < 0) -1 else 1

      dtr <- tr %>% select(all_of(c(o$ycol, preds))) %>% drop_na()
      dte <- te %>% select(all_of(c(o$ycol, preds))) %>% drop_na()
      if (nrow(dtr) < 50 || nrow(dte) < 50) next
      if (length(unique(dtr[[o$ycol]])) < 2 || length(unique(dte[[o$ycol]])) < 2) next
      for (pv in preds) { dtr[[pv]] <- sgn * dtr[[pv]]; dte[[pv]] <- sgn * dte[[pv]] }

      f <- as.formula(paste(o$ycol, "~", paste(preds, collapse = " + ")))
      m <- tryCatch(suppressWarnings(glm(f, data = dtr, family = binomial)),
                    error = function(e) NULL)
      if (is.null(m)) next

      p_tr <- suppressWarnings(predict(m, newdata = dtr, type = "response"))
      cut  <- youden_cut(dtr[[o$ycol]], p_tr)
      if (is.na(cut)) next
      p_te <- suppressWarnings(predict(m, newdata = dte, type = "response"))

      r_tr <- pROC::roc(dtr[[o$ycol]], p_tr, quiet = TRUE, direction = "<")
      r_te <- pROC::roc(dte[[o$ycol]], p_te, quiet = TRUE, direction = "<")
      auc_tr <- as.numeric(pROC::auc(r_tr))
      auc_te <- as.numeric(pROC::auc(r_te))
      ci <- as.numeric(pROC::ci.auc(r_te))
      dm <- diag_metrics(dte[[o$ycol]], as.numeric(p_te >= cut))
      cal <- calibration(dte[[o$ycol]], p_te)
      # 这里只描述训练 AUC 点估计与验证 AUC CI 的位置关系；不是正式的
      # 训练—验证 AUC 差检验，不能据此标注“显著”或“等效/保持”。
      ci_relation <- if (auc_tr > ci[3]) "训练AUC高于验证CI上界"
                     else if (auc_tr < ci[1]) "训练AUC低于验证CI下界"
                     else "训练AUC位于验证CI内"

      cat(sprintf("  [%s → %s] %-9s %-14s 阈值=%.3f\n", tr_name, te_name, o$outcome, mm$label, cut))
      cat(sprintf("        训练 AUC=%.3f（n=%d，患病率%.1f%%） → 验证 AUC=%.3f [%.3f, %.3f]（n=%d，患病率%.1f%%）｜ΔAUC=%+.3f｜%s\n",
                  auc_tr, nrow(dtr), mean(dtr[[o$ycol]]) * 100,
                  auc_te, ci[1], ci[3], nrow(dte), mean(dte[[o$ycol]]) * 100,
                  auc_te - auc_tr, ci_relation))
      cat(sprintf("        Sens=%.3f Spec=%.3f PPV=%.3f NPV=%.3f LR+=%.2f LR−=%.2f 准确率=%.3f\n",
                  dm["sens"], dm["spec"], dm["ppv"], dm["npv"], dm["lr_pos"], dm["lr_neg"], dm["acc"]))
      cat(sprintf("        校准斜率=%.3f（理想 1）｜校准截距=%.3f（理想 0）｜HL p=%.3f（大样本下敏感，仅参考）\n",
                  cal["cal_slope"], cal["cal_intercept"], cal["hl_p"]))

      rows[[length(rows) + 1]] <- data.frame(
        train_cohort = tr_name, test_cohort = te_name,
        outcome = o$outcome, model = mm$label,
        n_train = nrow(dtr), n_test = nrow(dte),
        prev_train = round(mean(dtr[[o$ycol]]), 4), prev_test = round(mean(dte[[o$ycol]]), 4),
        cutpoint = round(cut, 4),
        auc_train = round(auc_tr, 4),
        auc_test = round(auc_te, 4), auc_test_lo = round(ci[1], 4), auc_test_hi = round(ci[3], 4),
        delta_auc = round(auc_te - auc_tr, 4),
        auc_ci_relation = ci_relation,
        sens = round(dm["sens"], 4), spec = round(dm["spec"], 4),
        ppv = round(dm["ppv"], 4), npv = round(dm["npv"], 4),
        lr_pos = round(dm["lr_pos"], 3), lr_neg = round(dm["lr_neg"], 3),
        acc = round(dm["acc"], 4),
        cal_slope = round(cal["cal_slope"], 3),
        cal_intercept = round(cal["cal_intercept"], 3),
        hl_p = fmt_p(cal["hl_p"]),
        stringsAsFactors = FALSE)
    }
  }
}

if (length(rows)) {
  tab <- do.call(rbind, rows)
  save_tab(tab, "T20_external_validation")

  cat("\n  ★ 汇总（四 α 联合模型的跨队列表现）\n")
  sub <- tab[tab$model == "四 α 联合", ]
  if (nrow(sub)) print(as.data.frame(sub[, c("train_cohort","test_cohort","outcome",
                                             "auc_train","auc_test","delta_auc","auc_ci_relation")]),
                       row.names = FALSE)
  cat("\n  ★ 判读要点：\n")
  cat("    · 训练与验证 AUC 及其区间需并列报告；点估计落入对方 CI 不构成差异检验或等效证据\n")
  cat("    · PPV/NPV 会随患病率与切点变化，部署前需在目标场景重新评价，而非直接沿用本样本阈值\n")
  cat("    · 校准斜率与截距描述验证样本中的校准偏离，不能单独归因于某一种人群机制\n")
  cat("    · 单变量与多变量均先概率化；验证 AUC、CI 与校准使用同一连续预测值\n")
} else {
  cat("  ⚠ 未生成任何行（数据不足）\n")
}
