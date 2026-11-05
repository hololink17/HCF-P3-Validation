# =============================================================
# 10_analysis_25_replication_extension.R —— R11 复现核定与扩展（2026-10-02 新增）
# -------------------------------------------------------------
# 背景：2026-06-11 曾向另一篇论文（Front. Psychol. 17:1870659）的审稿人提交
#   私有补充材料（OSF 3k62h），内含 Prolific 初步结果。现役管线此前未覆盖其中若干项。
#   本模块把可复现/可补做的项**正式纳入管线**，逐条给出实测结果。
#   （判定与写作口径见 写论文/3k62h十条主张处置与偏离三重基准_v1_20261002.docx）
#
# 本模块保留公开 H_trend = α1² 的输出级内部一致性，以及通用下游统计：
# 肥胖 logistic、BMI 回归、数值锚点信度、四 α 潜剖面和外部分布比较。
# H_trend 仅为公开下游变换；定义式恒等不是结果，主检验为 α1² 对后台 H 的关联。
# 公开 L2 不含原始题目作答；相关分支在输入缺失时显式跳过。
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

set.seed(20261002)
suppressPackageStartupMessages(library(mclust))

ALPHAS <- c("alpha1", "alpha2", "alpha3", "alpha4")
ANCHORS <- paste0("anchor", 1:4)
COVARS <- c("age_group", "gender", "education")

# 两队列共有的 7 项自报疾病（见数无睡眠障碍阳性变量，故不纳入睡眠；心脏病补位）
DIS7 <- c("ibs_positive", "chronic_pain", "heart_disease", "hypertension",
          "diabetes", "hyperlipidemia", "obese")
DIS7_LAB <- c("IBS", "慢性痛", "心脏病", "高血压", "糖尿病", "高脂血症", "肥胖")

# ---------- 单队列还是双队列 ----------
L <- if (P3_COHORT %in% c("prolific", "jianshu")) load_cohorts(P3_COHORT) else load_both()
for (nm in names(L)) L[[nm]]$dep <- as.numeric(L[[nm]]$phq9_total >= 10)

# ---------- 共用工具 ----------
cronbach_alpha <- function(df) {
  df <- df[stats::complete.cases(df), , drop = FALSE]
  k <- ncol(df)
  if (nrow(df) < 10 || k < 2) return(c(alpha = NA_real_, n = nrow(df)))
  c(alpha = (k / (k - 1)) * (1 - sum(apply(df, 2, var)) / var(rowSums(df))), n = nrow(df))
}

fit_or <- function(d, y, x, covars = NULL, unit, lab) {
  dd <- d %>% select(all_of(c(y, x, covars))) %>% drop_na()
  if (nrow(dd) < 40 || length(unique(dd[[y]])) < 2) return(NULL)
  f <- as.formula(paste(y, "~", paste(c(x, covars), collapse = " + ")))
  m <- suppressWarnings(glm(f, data = dd, family = binomial))
  b <- coef(m)[x]
  ci <- tryCatch(suppressMessages(confint(m))[x, ], error = function(e) c(NA_real_, NA_real_))
  data.frame(
    y = y, alpha = lab,
    model = if (is.null(covars)) "crude" else "adjusted",
    unit = unit,
    OR = round(exp(b * unit), 4),
    lo = round(exp(ci[1] * unit), 4), hi = round(exp(ci[2] * unit), 4),
    p = format.pval(summary(m)$coefficients[x, 4], digits = 3),
    n = nrow(dd), event_rate = round(mean(dd[[y]]), 4),
    stringsAsFactors = FALSE)
}

# =============================================================
# 25.0　H_trend 输出级内部一致性
# =============================================================
cat("\n=== 25.0　H_trend = α1² 与后台 H（横断面内部一致性）===\n")
hrows <- list()
# 与 Python statsmodels 对齐：AIC 只计回归系数，不另计残差方差参数。
aic_ols <- function(m) -2 * as.numeric(logLik(m)) + 2 * length(coef(m))
for (nm in names(L)) {
  d <- L[[nm]]
  k <- is.finite(d$alpha1) & is.finite(d$H_score)
  a1sq <- (d$alpha1^2)[k]
  Hs <- d$H_score[k]
  r_rep <- cor(a1sq, Hs)
  mq <- lm(Hs ~ a1sq)
  ml <- lm(Hs ~ d$alpha1[k])
  cat(sprintf("  %s：cor([α1]², H) = %.4f；R² = %.4f；n = %d\n",
              nm, r_rep, r_rep^2, sum(k)))
  cat(sprintf("        平方模型 AIC = %.2f / 线性模型 AIC = %.2f\n",
              aic_ols(mq), aic_ols(ml)))
  hrows[[length(hrows) + 1]] <- data.frame(
    cohort = nm, n = sum(k),
    r_alpha1sq_vs_H = round(r_rep, 4), R2_alpha1sq_vs_H = round(r_rep^2, 4),
    AIC_quadratic = round(aic_ols(mq), 2), AIC_linear = round(aic_ols(ml), 2),
    note = "H_trend is the public downstream transformation alpha1^2; identity itself is not evidence",
    stringsAsFactors = FALSE)
}
save_tab(do.call(rbind, hrows), "T25_htrend_replication")

# =============================================================
# 25.1　α → 肥胖（logistic）
# =============================================================
cat("\n=== 25.1　α → 肥胖（BMI≥30，logistic）===\n")
orows <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  cat(sprintf("  --- %s（n=%d，肥胖率=%.3f）---\n", nm, nrow(d), mean(d$obese, na.rm = TRUE)))
  for (a in ALPHAS) {
    sdv <- sd(d[[a]], na.rm = TRUE)
    for (u in list(c("per +0.10", 0.10), c("per +1SD", sdv), c("per +1.00", 1.00))) {
      for (cv in list(NULL, COVARS)) {
        r <- fit_or(d, "obese", a, cv, as.numeric(u[2]), a)
        if (!is.null(r)) { r$cohort <- nm; r$unit <- u[1]; orows[[length(orows) + 1]] <- r }
      }
    }
    r <- fit_or(d, "obese", a, COVARS, 0.10, a)
    if (!is.null(r))
      cat(sprintf("    %s（adjusted, per+0.10）：OR=%.4f (%.3f–%.3f)  p=%s  n=%d\n",
                  a, r$OR, r$lo, r$hi, r$p, r$n))
  }
}
ob <- do.call(rbind, orows)
ob <- ob[, c("cohort", "y", "alpha", "model", "unit", "OR", "lo", "hi", "p", "n", "event_rate")]
save_tab(ob, "T25_obesity_or")

# =============================================================
# 25.2　H12：α → BMI（线性回归）
# =============================================================
cat("\n=== 25.2　H12：α → BMI ===\n")
brows <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  for (a in ALPHAS) {
    dd <- d %>% select(all_of(c("bmi", a, COVARS))) %>% drop_na()
    if (nrow(dd) < 40) next
    m0 <- lm(as.formula(paste("bmi ~", paste(COVARS, collapse = " + "))), data = dd)
    m1 <- lm(as.formula(paste("bmi ~", a, "+", paste(COVARS, collapse = " + "))), data = dd)
    av <- anova(m0, m1); cf <- summary(m1)$coefficients
    cat(sprintf("  %s %s：β=%.3f  p=%s  ΔR²=%+.4f (p=%s)  n=%d\n",
                nm, a, cf[a, 1], format.pval(cf[a, 4], digits = 3),
                summary(m1)$r.squared - summary(m0)$r.squared,
                format.pval(av$`Pr(>F)`[2], digits = 3), nrow(dd)))
    brows[[length(brows) + 1]] <- data.frame(
      cohort = nm, alpha = a, beta = round(cf[a, 1], 4), se = round(cf[a, 2], 4),
      t = round(cf[a, 3], 3), p = format.pval(cf[a, 4], digits = 3),
      R2_covars_only = round(summary(m0)$r.squared, 4),
      R2_full = round(summary(m1)$r.squared, 4),
      delta_R2 = round(summary(m1)$r.squared - summary(m0)$r.squared, 4),
      delta_R2_p = format.pval(av$`Pr(>F)`[2], digits = 3),
      n = nrow(dd), stringsAsFactors = FALSE)
  }
}
save_tab(do.call(rbind, brows), "T25_bmi_regression")

# =============================================================
# 25.3　锚题内部一致性（Cronbach α）多口径 + 收敛效度
# =============================================================
cat("\n=== 25.3　锚题内部一致性（Cronbach α）多口径 ===\n")
arows <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  subsets <- list(全队列 = rep(TRUE, nrow(d)))
  if ("data_source" %in% names(d) && any(d$data_source == "full", na.rm = TRUE))
    subsets[["完整版子集"]] <- !is.na(d$data_source) & d$data_source == "full"
  if (!all(ANCHORS %in% names(d)) ||
      !any(complete.cases(d[, ANCHORS, drop = FALSE]))) {
    cat(sprintf("  （%s 公开基底不含四个锚点题，跳过 25.3；不生成占位结果表）\n", nm))
    next
  }
  for (sn in names(subsets)) {
    idx <- subsets[[sn]]
    sub <- as.data.frame(lapply(d[idx, ANCHORS, drop = FALSE],
                                function(z) suppressWarnings(as.numeric(z))))
    ca <- cronbach_alpha(sub)
    rs <- vapply(1:4, function(i) {
      x <- suppressWarnings(as.numeric(d[[ALPHAS[i]]]))[idx]
      y <- suppressWarnings(as.numeric(d[[ANCHORS[i]]]))[idx]
      k <- is.finite(x) & is.finite(y)
      if (sum(k) < 10) NA_real_ else round(cor(x[k], y[k]), 3)
    }, numeric(1))
    cat(sprintf("  %s / %s：n=%d  Cronbach α=%.4f ｜ 收敛 r：α₁=%.3f α₂=%.3f α₃=%.3f α₄=%.3f\n",
                nm, sn, ca["n"], ca["alpha"], rs[1], rs[2], rs[3], rs[4]))
    arows[[length(arows) + 1]] <- data.frame(
      cohort = nm, subset = sn, n = unname(ca["n"]),
      cronbach_alpha = round(unname(ca["alpha"]), 4),
      r_alpha1_anchor = rs[1], r_alpha2_anchor = rs[2],
      r_alpha3_anchor = rs[3], r_alpha4_anchor = rs[4],
      stringsAsFactors = FALSE)
  }
}
if (length(arows)) save_tab(do.call(rbind, arows), "T25_anchor_reliability")
cat("  ★ 判读纪律：\n")
cat("    · 报告锚题信度时**四个口径值必须同现**并标明子样本，禁只报最低值\n")
cat("    · 4 题分测 α₁–α₄ 四个不同构念，Cronbach α 前提（同质题目）不成立\n")
cat("      ⇒ α 低属**结构决定**，不得据此推断「用户不能可靠自评」\n")

# =============================================================
# 25.4　LPA 四 α 潜剖面（mclust）
# =============================================================
cat("\n=== 25.4　LPA 四 α 潜剖面（mclust, modelNames=VVI）===\n")
frows <- list(); prows <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  need <- c(ALPHAS, DIS7, "H_score", "phq9_total")
  dd <- d %>% select(all_of(need)) %>% drop_na()
  if (nrow(dd) < 100) { cat(sprintf("  （%s 样本不足，跳过）\n", nm)); next }
  X <- scale(as.matrix(dd[, ALPHAS]))
  # 逐 G 拟合（避开 mclustBIC 返回矩阵在不同版本下方向不一致的坑）
  fits <- list()
  for (g in 2:6) fits[[as.character(g)]] <- tryCatch(Mclust(X, G = g, modelNames = "VVI"),
                                                     error = function(e) NULL)
  bv <- vapply(names(fits), function(k) {
    f <- fits[[k]]; if (is.null(f)) NA_real_ else as.numeric(f$bic)
  }, numeric(1))
  names(bv) <- names(fits)
  if (all(is.na(bv))) { cat(sprintf("  （%s LPA 拟合全部失败，跳过）\n", nm)); next }
  best_G <- as.integer(names(bv)[which.max(bv)])
  cat(sprintf("  --- %s（n=%d）---\n", nm, nrow(dd)))
  for (g in names(bv))
    cat(sprintf("    G=%s  BIC=%.1f  %s\n", g, bv[g], ifelse(as.integer(g) == best_G, "← BIC 最优", "")))
  fit <- fits[[as.character(best_G)]]
  Z <- fit$z
  ent <- 1 - sum(-Z * log(pmax(Z, 1e-12))) / (nrow(Z) * log(fit$G))
  dd$cls <- fit$classification
  dd$disease_count <- rowSums(sapply(DIS7, function(x) as.numeric(dd[[x]] == 1)))
  cat(sprintf("    BIC 最优类数 = %d；分类熵 = %.3f\n", best_G, ent))
  for (g in sort(unique(dd$cls))) {
    s <- dd[dd$cls == g, ]
    cat(sprintf("    类%d：n=%3d（%4.1f%%）  H=%.3f  PHQ-9=%.2f  病数=%.2f  α=(%.3f,%.3f,%.3f,%.3f)\n",
                g, nrow(s), 100 * nrow(s) / nrow(dd), mean(s$H_score),
                mean(s$phq9_total), mean(s$disease_count),
                mean(s$alpha1), mean(s$alpha2), mean(s$alpha3), mean(s$alpha4)))
    prows[[length(prows) + 1]] <- data.frame(
      cohort = nm, best_G = best_G, class = g, n = nrow(s),
      pct = round(100 * nrow(s) / nrow(dd), 2),
      alpha1 = round(mean(s$alpha1), 4), alpha2 = round(mean(s$alpha2), 4),
      alpha3 = round(mean(s$alpha3), 4), alpha4 = round(mean(s$alpha4), 4),
      disease_count = round(mean(s$disease_count), 3),
      phq9 = round(mean(s$phq9_total), 3), H = round(mean(s$H_score), 4),
      stringsAsFactors = FALSE)
  }
  f4 <- fits[["4"]]
  if (!is.null(f4)) {
    c4 <- sort(table(f4$classification), decreasing = TRUE)
    cat(sprintf("    [强制 G=4] 各类 n = %s（BIC=%.1f）\n",
                paste(c4, collapse = " / "), bv["4"]))
    frows[[length(frows) + 1]] <- data.frame(
      cohort = nm, forced_G4_BIC = round(as.numeric(bv["4"]), 1),
      forced_G4_sizes = paste(c4, collapse = "/"),
      best_G = best_G, best_BIC = round(max(bv, na.rm = TRUE), 1),
      entropy_best = round(ent, 4), n = nrow(dd), stringsAsFactors = FALSE)
  }
}
if (length(frows)) save_tab(do.call(rbind, frows), "T25_lpa_fit")
if (length(prows)) save_tab(do.call(rbind, prows), "T25_lpa_profiles")

# =============================================================
# 25.5　跨疾病 top-1 选项一致性
# 该分支仅在受控题目级输入存在时运行；公开 L2 显式跳过。
# =============================================================
cat("\n=== 25.5　跨疾病 top-1 选项一致性（仅用原始作答列）===\n")
ITEM_COLS <- paste0("hcf_", 5:24)
# 每个队列各自的原始作答来源（按队列取，不按运行模式取）
RAW_SRC <- list(
  prolific = list(f = base_path("p3_base_gold"),
                  pid = "participant_token"),
  jianshu  = list(f = base_path("p3_base_jianshu"),
                  pid = "作答ID")
)
trows <- list()
for (nm in names(L)) {
  src <- RAW_SRC[[nm]]
  if (is.null(src) || !file.exists(src$f)) { cat(sprintf("  （%s 原始作答文件缺失，跳过）\n", nm)); next }
  raw <- readr::read_csv(src$f, col_types = readr::cols(.default = readr::col_character()),
                         show_col_types = FALSE)
  id_col <- intersect(c(src$pid, "participant_token", "pid", "作答ID"), names(raw))[1]
  if (is.na(id_col) || !any(ITEM_COLS %in% names(raw))) {
    cat(sprintf("  （%s 公开基底不含原始作答或可用 ID，跳过 25.5）\n", nm))
    next
  }
  have <- c(id_col, intersect(ITEM_COLS, names(raw)))
  raw <- raw %>% select(all_of(have)) %>% distinct(.data[[id_col]], .keep_all = TRUE)
  d <- L[[nm]]
  m <- d %>% select(pid, all_of(DIS7), dep) %>% distinct(pid, .keep_all = TRUE)
  j <- dplyr::inner_join(m, raw, by = c("pid" = id_col))
  if (nrow(j) < 50) { cat(sprintf("  （%s 与原始作答对齐不足，跳过 n=%d）\n", nm, nrow(j))); next }
  ds_cols <- c(DIS7, "dep"); ds_lab <- c(DIS7_LAB, "抑郁")
  agree <- c()
  for (it in intersect(ITEM_COLS, names(j))) {
    modes <- vapply(ds_cols, function(ds) {
      v <- j[[it]][!is.na(j[[ds]]) & j[[ds]] == 1]
      v <- v[!is.na(v) & v != ""]
      if (!length(v)) return(NA_character_)
      names(sort(table(v), decreasing = TRUE))[1]
    }, character(1))
    modes <- modes[!is.na(modes)]
    if (!length(modes)) next
    top <- names(sort(table(modes), decreasing = TRUE))[1]
    agree[it] <- sum(modes == top)
    trows[[length(trows) + 1]] <- data.frame(
      cohort = nm, item = it, modal_option = substr(top, 1, 40),
      agree_diseases = sum(modes == top), total_diseases = length(modes),
      agree_pct = round(100 * sum(modes == top) / length(modes), 1),
      stringsAsFactors = FALSE)
  }
  if (length(agree))
    cat(sprintf("  %s：%d 题中 %d 题 top-1 选项在全部 %d 病种间完全一致；平均一致率 %.1f%%\n",
                nm, length(agree), sum(agree == length(ds_cols)), length(ds_cols),
                100 * mean(agree / length(ds_cols))))
}
if (length(trows)) save_tab(do.call(rbind, trows), "T25_top1_option")

cat("\n=== 25 号模块完成 ===\n")
cat("  ★ 本模块覆盖 3k62h 十条主张中的第 ②③⑤⑥⑦⑧⑨ 条；\n")
cat("    处置与写作口径见：写论文/3k62h十条主张处置与偏离三重基准_v1_20261002.docx\n")
