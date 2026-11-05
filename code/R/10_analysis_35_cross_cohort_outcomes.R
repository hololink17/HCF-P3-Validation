# =============================================================
# 10_analysis_35_cross_cohort_outcomes.R —— 跨队列结局对照（补齐）
# 用户要求：肥胖 / 慢性疼痛 / IBS(切点·剂量反应) / NNT / U型 / 决策树
#   的见数↔Prolific 跨队列对比。
# 与模块 11 的分工：模块 11 只比 α→IBS 预测 AUC + 分型/患病率；
#   本模块补齐「具体临床结局」的跨队列对照，全部用 load_both() 强制双队列、
#   从原始数据统一重算，保证两队列口径完全一致（不依赖单队列产物命名差异）。
# 产出（写入当前 OUT，jianshu/prolific 两次运行各写一份、内容相同）：
#   T35cc_obesity_or.csv         肥胖率 + α→肥胖 crude/adjusted OR（per +0.10）
#   T35cc_chronicpain_or.csv     慢性痛率 + α→慢性痛 OR
#   T35cc_ibs_dose_response.csv   α₃ 三分位 IBS 率 + 趋势OR / 连续OR
#   T35cc_nnt.csv                 IBS 筛查 NNT（α₃ 逻辑回归，阈值 0.15–0.40）
#   T35cc_ibs_youden_cutoff.csv   IBS 最优 Youden 切点（α₃ 与 H_score）
#   T35cc_chronicpain_ushape.csv 慢性痛 ~ α₃ 二次型（U型）检验
#   T35cc_decision_tree_ibs.csv   IBS 决策树变量重要性（跨队列谁先分裂）
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_both()   # 强制读入两队列
ALPHAS <- c("alpha1", "alpha2", "alpha3", "alpha4")
COVARS <- c("age_group", "gender", "education")
rows <- list(); cat("\n=== 模块35 跨队列结局对照（load_both 强制双队列）===\n")

# 统一派生 obese（BMI≥30），与单队列模块25一致
for (nm in names(L)) {
  d <- L[[nm]]
  d$obese <- as.numeric(d$bmi >= 30)
  L[[nm]] <- d
}
cov_ok <- function(d) all(COVARS %in% names(d)) && all(sapply(d[COVARS], function(x) length(unique(na.omit(x))) > 1))

# ---------- 通用 OR 拟合（与模块25 fit_or 同口径）----------
fit_or <- function(d, y, x, covars = NULL, unit) {
  dd <- d %>% select(all_of(c(y, x, covars))) %>% drop_na()
  if (nrow(dd) < 40 || length(unique(dd[[y]])) < 2) return(NULL)
  f <- as.formula(paste(y, "~", paste(c(x, covars), collapse = " + ")))
  m <- suppressWarnings(glm(f, data = dd, family = binomial))
  b <- coef(m)[x]
  ci <- tryCatch(suppressMessages(confint(m))[x, ], error = function(e) c(NA_real_, NA_real_))
  data.frame(y = y, alpha = x,
             model = if (is.null(covars)) "crude" else "adjusted",
             unit = sprintf("per +%.2f", unit),
             OR = round(exp(b * unit), 4),
             lo = round(exp(ci[1] * unit), 4), hi = round(exp(ci[2] * unit), 4),
             p = format.pval(summary(m)$coefficients[x, 4], digits = 3),
             n = nrow(dd), event_rate = round(mean(dd[[y]]), 4),
             stringsAsFactors = FALSE)
}

# =============================================================
# 1) 肥胖 & 2) 慢性痛 OR（跨队列）
# =============================================================
or_tbl <- function(y, lab) {
  out <- list()
  for (nm in names(L)) {
    d <- L[[nm]]
    er <- mean(d[[y]], na.rm = TRUE)
    cat(sprintf("  [%s] %s 率=%.3f (n=%d)\n", nm, lab, er, nrow(d)))
    for (a in ALPHAS) {
      r0 <- fit_or(d, y, a, NULL, 0.10)          # crude
      if (!is.null(r0)) { r0$cohort <- nm; out[[length(out)+1]] <- r0 }
      if (cov_ok(d)) {
        r1 <- fit_or(d, y, a, COVARS, 0.10)       # adjusted
        if (!is.null(r1)) { r1$cohort <- nm; out[[length(out)+1]] <- r1 }
      }
    }
  }
  do.call(rbind, out)
}
ob_tbl <- or_tbl("obese", "肥胖"); ob_tbl <- ob_tbl[, c("cohort","y","alpha","model","unit","OR","lo","hi","p","n","event_rate")]
save_tab(ob_tbl, "T35cc_obesity_or")
cp_tbl <- or_tbl("chronic_pain", "慢性痛"); cp_tbl <- cp_tbl[, c("cohort","y","alpha","model","unit","OR","lo","hi","p","n","event_rate")]
save_tab(cp_tbl, "T35cc_chronicpain_or")
cat("  已写 T35cc_obesity_or / T35cc_chronicpain_or\n")

# =============================================================
# 3) IBS 剂量反应（α₃ 三分位 + 趋势OR / 连续OR）
# =============================================================
dr <- list()
for (nm in names(L)) {
  d <- L[[nm]] %>% filter(is.finite(alpha3), is.finite(ibs_positive))
  q <- ntile(d$alpha3, 3)
  d$q3 <- q
  rates <- tapply(d$ibs_positive, d$q3, mean, na.rm = TRUE)
  ns    <- tapply(d$ibs_positive, d$q3, length)
  # 趋势 OR（每分位）
  dd <- d; dd$score <- as.numeric(dd$q3) - 1
  mt <- suppressWarnings(glm(ibs_positive ~ score, data = dd, family = binomial))
  trend_or <- round(exp(coef(mt)["score"]), 3); trend_p <- format.pval(summary(mt)$coefficients["score","Pr(>|z|)"], digits = 3)
  # 连续 OR（每 +0.10 α₃）
  mc <- suppressWarnings(glm(ibs_positive ~ alpha3, data = dd, family = binomial))
  cont_or <- round(exp(coef(mc)["alpha3"] * 0.10), 3); cont_p <- format.pval(summary(mc)$coefficients["alpha3","Pr(>|z|)"], digits = 3)
  # Q1 vs Q3 率差
  qdiff <- round(100*(rates[3] - rates[1]), 1)
  cat(sprintf("  [%s] α₃ 三分位 IBS率 Q1=%.1f%% Q2=%.1f%% Q3=%.1f%%  趋势OR/分位=%.3f(p=%s) 连续OR/+0.10=%.3f(p=%s)\n",
              nm, 100*rates[1], 100*rates[2], 100*rates[3], trend_or, trend_p, cont_or, cont_p))
  dr[[length(dr)+1]] <- data.frame(cohort=nm, analysis="α₃三分位IBS率",
    Q1_rate_pct=round(100*rates[1],1), Q1_n=ns[1],
    Q2_rate_pct=round(100*rates[2],1), Q2_n=ns[2],
    Q3_rate_pct=round(100*rates[3],1), Q3_n=ns[3],
    trend_OR_per_tertile=trend_or, trend_p=trend_p,
    cont_OR_per_0.10=cont_or, cont_p=cont_p,
    Q3_minus_Q1_rate_pct=qdiff, stringsAsFactors=FALSE)
}
save_tab(do.call(rbind, dr), "T35cc_ibs_dose_response")
cat("  已写 T35cc_ibs_dose_response\n")

# =============================================================
# 4) NNT（IBS 筛查，α₃ 逻辑回归，阈值 0.15–0.40，与模块06同口径）
# =============================================================
nnt_at <- function(y, p, pt) {
  sel <- p >= pt
  if (sum(sel) == 0) return(c(nnt=NA, n_sel=0, arr=NA))
  arr <- mean(y[sel]) - mean(y[!sel])
  c(nnt = 1/abs(arr), n_sel = sum(sel), arr = arr)
}
nnt_rows <- list()
for (nm in names(L)) {
  d <- L[[nm]] %>% filter(is.finite(alpha3), is.finite(ibs_positive))
  m <- suppressWarnings(glm(ibs_positive ~ alpha3, data = d, family = binomial))
  p <- predict(m, type = "response")
  r <- pROC::roc(d$ibs_positive, p, quiet = TRUE)
  pt_best <- as.numeric(pROC::coords(r, "best", ret = "threshold")[1,1])
  a3_at <- (log(pt_best/(1-pt_best)) - coef(m)[1]) / coef(m)[2]
  cat(sprintf("  [%s] IBS 筛查 Youden切点 p=%.4f → α₃=%.4f  NNT@0.15=?", nm, pt_best, a3_at))
  for (t in c(0.15, 0.20, 0.25, 0.30, 0.40)) {
    v <- nnt_at(d$ibs_positive, p, t)
    cat(sprintf(" %.1f", ifelse(is.na(v["nnt"]), NA, v["nnt"])))
    nnt_rows[[length(nnt_rows)+1]] <- data.frame(cohort=nm, threshold=t,
      NNT=round(v["nnt"],2), n_selected=v["n_sel"], ARR=round(v["arr"],4),
      youden_pt=round(pt_best,4), a3_at_youden=round(a3_at,4), stringsAsFactors=FALSE)
  }
  cat("\n")
}
save_tab(do.call(rbind, nnt_rows), "T35cc_nnt")
cat("  已写 T35cc_nnt\n")

# =============================================================
# 5) IBS 最优 Youden 切点（α₃ 与 H_score）
# =============================================================
cut_rows <- list()
for (pred in c("alpha3","H_score")) {
  for (nm in names(L)) {
    d <- L[[nm]] %>% filter(is.finite(.data[[pred]]), is.finite(ibs_positive))
    r <- pROC::roc(d$ibs_positive, d[[pred]], quiet = TRUE)
    best <- pROC::coords(r, "best", ret = c("threshold","sensitivity","specificity"))
    cut_rows[[length(cut_rows)+1]] <- data.frame(cohort=nm, predictor=pred, outcome="ibs_positive",
      cutoff=round(as.numeric(best["threshold"]),4),
      sens=round(as.numeric(best["sensitivity"]),4),
      spec=round(as.numeric(best["specificity"]),4),
      youden=round(as.numeric(best["sensitivity"])+as.numeric(best["specificity"])-1,4),
      stringsAsFactors=FALSE)
    cat(sprintf("  [%s] %s→IBS 切点=%.4f (sens=%.2f spec=%.2f)\n", nm, pred,
                as.numeric(best["threshold"]), as.numeric(best["sensitivity"]), as.numeric(best["specificity"])))
  }
}
save_tab(do.call(rbind, cut_rows), "T35cc_ibs_youden_cutoff")
cat("  已写 T35cc_ibs_youden_cutoff\n")

# =============================================================
# 6) 慢性痛 U型（α₃ 二次型）
# =============================================================
ush <- list()
for (nm in names(L)) {
  d <- L[[nm]] %>% filter(is.finite(alpha3), is.finite(chronic_pain))
  m0 <- suppressWarnings(glm(chronic_pain ~ alpha3, data = d, family = binomial))
  m1 <- suppressWarnings(glm(chronic_pain ~ alpha3 + I(alpha3^2), data = d, family = binomial))
  q <- ntile(d$alpha3, 4); d$q4 <- q
  qr <- tapply(d$chronic_pain, d$q4, mean, na.rm = TRUE)
  lin_p <- format.pval(summary(m0)$coefficients["alpha3","Pr(>|z|)"], digits = 3)
  quad_p <- format.pval(summary(m1)$coefficients["I(alpha3^2)","Pr(>|z|)"], digits = 3)
  shape <- if (as.numeric(quad_p) < 0.05) "存在显著非线性(U型/J型)" else "线性/单调(无显著U型)"
  cat(sprintf("  [%s] 慢性痛~α₃ 线性p=%s 二次型p=%s → %s  (Q1–Q4率: %.0f%%→%.0f%%)\n",
              nm, lin_p, quad_p, shape, 100*qr[1], 100*qr[4]))
  ush[[length(ush)+1]] <- data.frame(cohort=nm, disease="慢性痛",
    linear_p=lin_p, quadratic_p=quad_p, shape=shape,
    Q1_rate_pct=round(100*qr[1],1), Q2_rate_pct=round(100*qr[2],1),
    Q3_rate_pct=round(100*qr[3],1), Q4_rate_pct=round(100*qr[4],1), stringsAsFactors=FALSE)
}
save_tab(do.call(rbind, ush), "T35cc_chronicpain_ushape")
cat("  已写 T35cc_chronicpain_ushape\n")

# =============================================================
# 7) IBS 决策树变量重要性（跨队列谁先分裂）
# =============================================================
tr <- list()
for (nm in names(L)) {
  d <- L[[nm]] %>% select(all_of(c("ibs_positive", ALPHAS))) %>% drop_na()
  fit <- suppressWarnings(rpart::rpart(as.factor(ibs_positive) ~ ., data = d,
    control = rpart::rpart.control(cp = 0.02, minsplit = 20)))
  imp <- sort(fit$variable.importance, decreasing = TRUE)
  if (is.null(imp) || length(imp) == 0) {
    top <- NA_character_; rank <- "无有效分裂（树为空）"
    cat(sprintf("  [%s] IBS 决策树：无有效分裂（cp/minSplit 过严或预测力弱）\n", nm))
  } else {
    top <- names(imp)[1]
    rank <- paste(names(imp), round(as.numeric(imp), 2), collapse = " > ")
    cat(sprintf("  [%s] IBS 决策树 top变量=%s (重要性: %s)\n", nm, top, rank))
  }
  tr[[length(tr)+1]] <- data.frame(cohort=nm, top_split=top,
    importance_rank=rank, stringsAsFactors=FALSE)
}
save_tab(do.call(rbind, tr), "T35cc_decision_tree_ibs")
cat("  已写 T35cc_decision_tree_ibs\n")

cat("=== 模块35 跨队列结局对照完成 ===\n")
