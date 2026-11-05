# =============================================================
# 10_analysis_36_seven_factors.R —— 7 因子主分析（α₁/α₂/α₄ + a/b/c/d）
# -------------------------------------------------------------
# 口径声明：本轮分析对象 = α₁、α₂、α₄ + 载体 a、b、c、d；
# α₃ 的独立结构另见专门模块。所有读数均来自单一固定后台管线。
#
# 本模块产出：
#   ① T36_7factor_descriptives.csv —— 7 因子描述统计（两平台 + 合并）
#   ② T36_7factor_correlations.csv —— 7 因子相关矩阵（两平台）
#   ③ T36_7factor_outcomes.csv —— 7 因子与效标（PHQ-9/GAD-7/IBS/BMI/疾病）的关联
#   ④ T36_7factor_auc.csv —— 单因子 AUC（含 95%CI）
#   ⑤ T36_7factor_composite.csv —— 7 因子复合分的增量效度（IBS  Logistic）
#   ⑥ T36_7factor_alpha3_note.csv —— α₃ 单独说明（不在 7 因子里）
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

F7 <- c("alpha1", "alpha2", "alpha4", "carrier_a", "carrier_b", "carrier_c", "carrier_d")
F7_ZH <- c(alpha1 = "α₁ 觉察", alpha2 = "α₂ 接纳", alpha4 = "α₄ 行动",
           carrier_a = "carrier_a", carrier_b = "carrier_b",
           carrier_c = "carrier_c", carrier_d = "carrier_d")

cat("\n", strrep("=", 74), "\n", sep = "")
cat("  模块 36：7 因子分析（α₁/α₂/α₄ + a/b/c/d；单一固定管线）\n")
cat(strrep("=", 74), "\n", sep = "")

L <- load_both()          # 强制双队列（两平台同规则）
cohorts <- list(Prolific = L$prolific, 见数 = L$jianshu)

# ---------- ① 描述统计 ----------
desc_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  for (v in F7) {
    x <- d[[v]]
    x <- x[is.finite(x)]
    if (!length(x)) next
    desc_rows[[length(desc_rows) + 1]] <- data.frame(
      cohort = cn, factor = v, factor_zh = unname(F7_ZH[v]),
      n = length(x), mean = mean(x), sd = sd(x),
      median = median(x), q25 = quantile(x, .25), q75 = quantile(x, .75),
      min = min(x), max = max(x),
      stringsAsFactors = FALSE)
  }
}
desc <- do.call(rbind, desc_rows)
cat("\n=== 一、7 因子描述统计 ===\n")
print(desc[, c("cohort", "factor_zh", "n", "mean", "sd", "median", "min", "max")],
      row.names = FALSE)
save_tab(desc, "T36_7factor_descriptives")

# ---------- ② 相关矩阵 ----------
cor_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  m <- as.matrix(d[, F7])
  m <- m[complete.cases(m), , drop = FALSE]
  if (nrow(m) < 20) next
  cm <- suppressWarnings(cor(m, use = "pairwise.complete.obs"))
  for (i in seq_along(F7)) for (j in seq_along(F7)) {
    if (i >= j) next
    tt <- suppressWarnings(cor.test(m[, i], m[, j]))
    cor_rows[[length(cor_rows) + 1]] <- data.frame(
      cohort = cn, f1 = F7[i], f2 = F7[j],
      r = tt$estimate, p = tt$p.value, n = sum(complete.cases(m[, c(i, j)])),
      stringsAsFactors = FALSE)
  }
}
cors <- do.call(rbind, cor_rows)
cat("\n=== 二、7 因子两两相关 ===\n")
for (cn in unique(cors$cohort)) {
  cat(sprintf("\n--- %s ---\n", cn))
  print(cors[cors$cohort == cn, c("f1", "f2", "r", "p", "n")], row.names = FALSE)
}
save_tab(cors, "T36_7factor_correlations")

# ---------- ③ 7 因子与效标关联 ----------
out_rows <- list()
OUTCOMES <- list(
  phq9_total = c("PHQ-9 总分", "numeric"),
  gad7_total = c("GAD-7 总分", "numeric"),
  bmi        = c("BMI", "numeric"),
  ibs_positive = c("IBS 阳性", "binary"),
  depression = c("抑郁筛查阳性(PHQ≥10)", "binary"),
  chronic_pain = c("慢性痛", "binary"),
  hypertension = c("高血压", "binary"),
  diabetes = c("糖尿病", "binary"),
  hyperlipidemia = c("高血脂", "binary"),
  heart_disease = c("心脏病", "binary")
)
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  for (v in F7) {
    for (oc in names(OUTCOMES)) {
      if (!oc %in% names(d)) next
      x <- d[[v]]; y <- d[[oc]]
      ok <- is.finite(x) & !is.na(y)
      if (sum(ok) < 30) next
      x <- x[ok]; y <- y[ok]
      if (OUTCOMES[[oc]][2] == "numeric") {
        tt <- suppressWarnings(cor.test(x, y))
        r <- tt$estimate; p <- tt$p.value
        out_rows[[length(out_rows) + 1]] <- data.frame(
          cohort = cn, factor = v, factor_zh = unname(F7_ZH[v]),
          outcome = oc, outcome_zh = OUTCOMES[[oc]][1], type = "pearson",
          r = r, p = p, n = length(x), stringsAsFactors = FALSE)
      } else {
        yb <- as.numeric(y)
        if (!all(is.finite(yb)) || !all(yb %in% c(0, 1))) stop("二分类结局必须编码为0/1: ", oc)
        tt <- suppressWarnings(cor.test(x, yb))
        # 点二列相关
        r_pb <- suppressWarnings(cor(x, yb))
        out_rows[[length(out_rows) + 1]] <- data.frame(
          cohort = cn, factor = v, factor_zh = unname(F7_ZH[v]),
          outcome = oc, outcome_zh = OUTCOMES[[oc]][1], type = "point-biserial",
          r = r_pb, p = tt$p.value, n = length(x), stringsAsFactors = FALSE)
      }
    }
  }
}
outs <- do.call(rbind, out_rows)
cat("\n=== 三、7 因子与效标关联（|r| 排序，各队列前 12）===\n")
for (cn in unique(outs$cohort)) {
  s <- outs[outs$cohort == cn, ]
  s$absr <- abs(s$r)
  s <- s[order(-s$absr), ]
  cat(sprintf("\n--- %s ---\n", cn))
  print(head(s[, c("factor_zh", "outcome_zh", "r", "p", "n")], 12), row.names = FALSE)
}
save_tab(outs[, c("cohort", "factor", "factor_zh", "outcome", "outcome_zh",
                   "type", "r", "p", "n")], "T36_7factor_outcomes")

# ---------- ④ 单因子 AUC ----------
auc_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  for (oc in c("ibs_positive", "depression")) {
    if (!oc %in% names(d)) next
    for (v in F7) {
      a <- auc_ci(d[[oc]], d[[v]], orient = TRUE)
      auc_rows[[length(auc_rows) + 1]] <- data.frame(
        cohort = cn, outcome = oc, factor = v, factor_zh = unname(F7_ZH[v]),
        auc = a[["auc"]], ci_lo = a[["lo"]], ci_hi = a[["hi"]],
        n = a[["n"]], pos = a[["pos"]], stringsAsFactors = FALSE)
    }
  }
}
aucs7 <- do.call(rbind, auc_rows)
cat("\n=== 四、7 因子单因子 AUC（IBS / 抑郁筛查）===\n")
print(aucs7, row.names = FALSE)
save_tab(aucs7, "T36_7factor_auc")

# ---------- ⑤ 增量效度：IBS 的 Logistic（逐因子加入） ----------
inc_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]] %>% select(all_of(c("ibs_positive", F7))) %>% drop_na()
  if (nrow(d) < 60 || length(unique(d$ibs_positive)) < 2) next
  base_m <- suppressWarnings(glm(ibs_positive ~ 1, data = d, family = binomial))
  b0 <- as.numeric(pROC::auc(pROC::roc(d$ibs_positive,
                                       predict(base_m, newdata = d, type = "response"),
                                       quiet = TRUE)))
  prev <- b0
  current <- character()
  for (v in F7) {
    current <- c(current, v)
    f <- as.formula(paste("ibs_positive ~", paste(current, collapse = " + ")))
    m <- suppressWarnings(glm(f, data = d, family = binomial))
    pr <- suppressWarnings(predict(m, newdata = d, type = "response"))
    a <- as.numeric(pROC::auc(pROC::roc(d$ibs_positive, pr, quiet = TRUE)))
    inc_rows[[length(inc_rows) + 1]] <- data.frame(
      cohort = cn, step = match(v, F7),
      added = v, added_zh = unname(F7_ZH[v]),
      auc_in = a, delta_auc = a - prev, n = nrow(d),
      stringsAsFactors = FALSE)
    prev <- a
  }
  # 第7步已经就是全7因子联合模型，不再追加一行重复模型。
}
inc <- do.call(rbind, inc_rows)
cat("\n=== 五、IBS 增量效度（按F7顺序累积加入的样本内 AUC；非预注册顺序）===\n")
for (cn in unique(inc$cohort)) {
  cat(sprintf("\n--- %s ---\n", cn))
  print(inc[inc$cohort == cn, ], row.names = FALSE)
}
save_tab(inc, "T36_7factor_composite")

# ---------- ⑥ α₃ 说明（本轮不纳入 7 因子） ----------
note_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  a3 <- d$alpha3; a3 <- a3[is.finite(a3)]
  note_rows[[length(note_rows) + 1]] <- data.frame(
    cohort = cn, n = length(a3), mean = mean(a3), sd = sd(a3),
    note = "α₃ 为独立后台输出读数，本轮 7 因子分析不含 α₃；此处仅作分析范围说明",
    stringsAsFactors = FALSE)
}
note <- do.call(rbind, note_rows)
cat("\n=== 六、α₃ 口径说明（本轮 7 因子不含 α₃）===\n")
print(note, row.names = FALSE)
save_tab(note, "T36_7factor_alpha3_note")

# ---------- ⑦ 因子独立性诊断（方法学关键）----------
# ⚠ 多个后台输出读数共享同一问卷输入信号，因此可能呈结构相关。
#   本公开脚本不描述或重构任何后台生成规则。
#   → 7 个读数彼此高度相关时，判别力无法区分「哪个读数在起作用」。
#   本节把这一结构性事实量化，供正文据实交代（不可把它误读成「7 个维度各自独立有效」）。
diag_rows <- list()
for (cn in names(cohorts)) {
  d <- cohorts[[cn]]
  m <- as.matrix(d[, F7]); m <- m[complete.cases(m), , drop = FALSE]
  if (nrow(m) < 20) next
  # 平均两两相关
  cm <- suppressWarnings(cor(m, use = "pairwise.complete.obs"))
  pair_r <- cm[upper.tri(cm)]
  # 主成分：前 1 主成分解释方差比例
  pc <- prcomp(m, center = TRUE, scale. = FALSE)
  ev <- (pc$sdev^2) / sum(pc$sdev^2)
  # 每个因子与「其余 6 因子平均分」的 R
  gm <- rowMeans(m)
  for (v in F7) {
    j <- match(v, F7)
    others <- rowMeans(m[, -j, drop = FALSE])
    r_v <- suppressWarnings(cor(m[, j], others))
    diag_rows[[length(diag_rows) + 1]] <- data.frame(
      cohort = cn, factor = v, factor_zh = unname(F7_ZH[v]),
      r_with_rest_mean = r_v,
      mean_r_all_pairs = mean(pair_r),
      pc1_var_explained = ev[1],
      pc2_var_explained = ev[2],
      n = nrow(m), stringsAsFactors = FALSE)
  }
}
diag7 <- do.call(rbind, diag_rows)
cat("\n=== 七、因子独立性诊断（重要方法学事实）===\n")
cat("   说明：这些后台读数共享同一问卷输入，不能据此假定其相互独立。\n")
cat("         两两相关与联合模型结果仅描述公开读数的统计关系。\n")
cat("         不得将读数关联反推或表述为后台生成规则。\n\n")
for (cn in unique(diag7$cohort)) {
  cat(sprintf("--- %s ---\n", cn))
  print(diag7[diag7$cohort == cn, ], row.names = FALSE)
  cat("\n")
}
save_tab(diag7, "T36_7factor_independence_diagnostic")

cat("\n  → 模块 36 完成：7 因子（α₁/α₂/α₄ + a/b/c/d），单一固定管线\n")
