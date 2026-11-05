# =============================================================
# 10_analysis_23_clinical_operability.R —— R4 【零成本补强】临床可操作性与已知组效度
# -------------------------------------------------------------
# 为什么要做：
#   ① 医生真正关心的不是 AUC，而是「我门诊里 10 个患者，这套工具能筛出几个真的」——
#      即 PPV / NPV。既有模块只报了 AUC / NNT / ARR，**没有 PPV / NPV 随患病率表**。
#   ② 已知组效度（known-groups validity）散落在 T7 / table1 等处，未整合为标准表。
# 本模块产出：
#   T23_ppv_npv            —— 按本队列标定的切点，在不同患病率下的 PPV / NPV / 似然比
#   T23_known_groups       —— α 分型 × PHQ-9 五级 / GAD-7 四级 / 疾病数分层
# 说明：切点与患病率均为**队列特异**，两队列各算各的（Prolific 与见数患病率不同）。
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

L1 <- load_cohorts(); d <- L1[[1]]; coh <- d$cohort[1]
cat(sprintf("\n=== 临床可操作性与已知组效度（%s，n=%d）===\n", coh, nrow(d)))

# 疾病计数（本模块内派生，避免依赖外部文件）
d$disease_count <- rowSums(cbind(d$chronic_pain, d$heart_disease, d$hypertension,
                                 d$diabetes, d$hyperlipidemia), na.rm = TRUE)
d$phq9_grade <- cut(d$phq9_total, c(-1, 4, 9, 14, 19, 27),
                    labels = c("无/极轻(0-4)", "轻度(5-9)", "中度(10-14)",
                               "中重度(15-19)", "重度(20-27)"))
d$gad7_grade <- cut(d$gad7_total, c(-1, 4, 9, 14, 21),
                    labels = c("无/极轻(0-4)", "轻度(5-9)", "中度(10-14)", "重度(15-21)"))

# ---------- 一、PPV / NPV 随患病率表 ----------
# 切点来源：本队列 α₃ 对 IBS 阳性的 Youden 最优点（方向：α₃ 越低风险越高）
cut_at <- function(y, x, lower_is_risk = TRUE) {
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]
  xs <- if (lower_is_risk) -x else x
  r <- pROC::roc(y, xs, quiet = TRUE)
  co <- pROC::coords(r, "best", best.method = "youden",
                     ret = c("threshold", "sensitivity", "specificity"), transpose = FALSE)
  list(cut_raw = as.numeric(co$threshold[1]) * (if (lower_is_risk) -1 else 1),
       sens = as.numeric(co$sensitivity[1]), spec = as.numeric(co$specificity[1]),
       auc = as.numeric(pROC::auc(r)))
}

ppv_npv <- function(sens, spec, prev) {
  ppv <- sens * prev / (sens * prev + (1 - spec) * (1 - prev))
  npv <- spec * (1 - prev) / (spec * (1 - prev) + (1 - sens) * prev)
  lrp <- sens / (1 - spec); lrn <- (1 - sens) / spec
  c(prev = prev, ppv = ppv, npv = npv, lr_pos = lrp, lr_neg = lrn,
    nns = ifelse(ppv - prev > 0, 1 / (ppv - prev), NA_real_))
}

prow <- list()
# 结局 1：IBS 阳性（以 α₃ 为筛查指标）
cd <- d %>% select(alpha3, ibs_positive) %>% drop_na()
if (nrow(cd) > 50 && length(unique(cd$ibs_positive)) == 2) {
  ct <- cut_at(cd$ibs_positive, cd$alpha3, lower_is_risk = TRUE)
  cat(sprintf("\n[IBS 阳性] α₃ Youden 切点=%.4f｜Sens=%.3f｜Spec=%.3f｜AUC=%.3f｜本队列患病率=%.1f%%\n",
              ct$cut_raw, ct$sens, ct$spec, ct$auc, mean(cd$ibs_positive) * 100))
  for (pv in c(0.05, 0.10, 0.15, 0.20, 0.25, 0.30)) {
    m <- ppv_npv(ct$sens, ct$spec, pv)
    prow[[length(prow) + 1]] <- data.frame(
      cohort = coh, outcome = "IBS 阳性", marker = "α₃",
      cutpoint = round(ct$cut_raw, 4), sens = round(ct$sens, 3), spec = round(ct$spec, 3),
      prevalence = sprintf("%.0f%%", pv * 100), ppv = round(m["ppv"], 3),
      npv = round(m["npv"], 3), lr_pos = round(m["lr_pos"], 2), lr_neg = round(m["lr_neg"], 2),
      nns = ifelse(is.na(m["nns"]), NA_real_, round(m["nns"], 1)), stringsAsFactors = FALSE)
    cat(sprintf("    患病率 %2.0f%% → PPV=%.3f  NPV=%.3f  LR+=%.2f  需筛查人数=%.1f\n",
                pv * 100, m["ppv"], m["npv"], m["lr_pos"], m["nns"]))
  }
}
# 结局 2：PHQ-9≥10（以 α₂ 为筛查指标——分层回归中 α₂ 贡献最大）
d$phq9_pos <- as.numeric(d$phq9_total >= 10)
cd2 <- d %>% select(alpha2, phq9_pos) %>% drop_na()
if (nrow(cd2) > 50 && length(unique(cd2$phq9_pos)) == 2) {
  ct2 <- cut_at(cd2$phq9_pos, cd2$alpha2, lower_is_risk = TRUE)
  cat(sprintf("\n[PHQ-9≥10] α₂ Youden 切点=%.4f｜Sens=%.3f｜Spec=%.3f｜AUC=%.3f｜本队列阳性率=%.1f%%\n",
              ct2$cut_raw, ct2$sens, ct2$spec, ct2$auc, mean(cd2$phq9_pos) * 100))
  for (pv in c(0.10, 0.20, 0.30, 0.40)) {
    m <- ppv_npv(ct2$sens, ct2$spec, pv)
    prow[[length(prow) + 1]] <- data.frame(
      cohort = coh, outcome = "PHQ-9≥10", marker = "α₂",
      cutpoint = round(ct2$cut_raw, 4), sens = round(ct2$sens, 3), spec = round(ct2$spec, 3),
      prevalence = sprintf("%.0f%%", pv * 100), ppv = round(m["ppv"], 3),
      npv = round(m["npv"], 3), lr_pos = round(m["lr_pos"], 2), lr_neg = round(m["lr_neg"], 2),
      nns = ifelse(is.na(m["nns"]), NA_real_, round(m["nns"], 1)), stringsAsFactors = FALSE)
    cat(sprintf("    患病率 %2.0f%% → PPV=%.3f  NPV=%.3f  LR+=%.2f  需筛查人数=%.1f\n",
                pv * 100, m["ppv"], m["npv"], m["lr_pos"], m["nns"]))
  }
}
if (length(prow)) save_tab(do.call(rbind, prow), "T23_ppv_npv")

# ---------- 二、已知组效度整合表 ----------
cramers_v <- function(tb) {
  tb <- as.matrix(tb)
  n <- sum(tb)
  if (n == 0) return(NA_real_)
  chi <- suppressWarnings(chisq.test(tb, correct = FALSE))
  k <- min(nrow(tb), ncol(tb))
  sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
}

krows <- list()
add_kg <- function(factor_name, level_tab) {
  cat(sprintf("\n[已知组效度] α 分型 × %s\n", factor_name))
  # ⚠ 先剔除全零的行/列：某队列若在某等级上人数为 0（如见数 PHQ-9 无「重度」组），
  #   全零列会使卡方期望频数为 0 → 贡献 0/0 → Cramér's V 变 NaN。
  level_tab <- level_tab[rowSums(level_tab) > 0, colSums(level_tab) > 0, drop = FALSE]
  if (nrow(level_tab) < 2 || ncol(level_tab) < 2) {
    cat("  ⚠ 有效行列不足（某等级无观测），跳过该族层\n"); return(invisible(NA_real_))
  }
  print(level_tab)
  v <- cramers_v(level_tab)
  ch <- suppressWarnings(chisq.test(as.matrix(level_tab), correct = FALSE))
  cat(sprintf("  Cramér's V=%.3f｜χ²(%d)=%.2f，p=%s\n",
              v, ch$parameter, as.numeric(ch$statistic), format.pval(ch$p.value, digits = 3)))
  # 落盘为**长表**（列联表各族层列数不同，宽表无法 rbind）
  df <- as.data.frame.matrix(level_tab)
  long <- data.frame(
    alpha_group = rep(rownames(df), times = ncol(df)),
    level       = rep(colnames(df), each = nrow(df)),
    n           = as.vector(as.matrix(df)),
    stringsAsFactors = FALSE)
  long$cramers_v <- round(v, 3)
  long$chi_p <- format.pval(ch$p.value, digits = 3)
  krows[[length(krows) + 1]] <<- cbind(
    data.frame(cohort = coh, stratifier = factor_name, stringsAsFactors = FALSE), long)
  invisible(v)
}

if (!all(is.na(d$alpha_group))) {
  add_kg("PHQ-9 五级", table(d$alpha_group, d$phq9_grade, useNA = "no"))
  add_kg("GAD-7 四级", table(d$alpha_group, d$gad7_grade, useNA = "no"))
  dc_lab <- cut(d$disease_count, c(-1, 0, 1, 2, Inf),
                labels = c("0 种", "1 种", "2 种", "3 种及以上"))
  add_kg("自报疾病数分层", table(d$alpha_group, dc_lab, useNA = "no"))
}

# 各分型组的症状均值与效应量（耗竭型 vs 高意识型）
cat("\n[已知组效度] 各分型组的症状均值（组间比较）\n")
g <- d %>% group_by(alpha_group) %>%
  summarise(n = n(),
            phq9_mean = round(mean(phq9_total, na.rm = TRUE), 2),
            phq9_sd = round(sd(phq9_total, na.rm = TRUE), 2),
            gad7_mean = round(mean(gad7_total, na.rm = TRUE), 2),
            gad7_sd = round(sd(gad7_total, na.rm = TRUE), 2),
            disease_mean = round(mean(disease_count, na.rm = TRUE), 2),
            .groups = "drop")
print(as.data.frame(g), row.names = FALSE)
save_tab(g, "T23_group_means")

# 效应量：耗竭型 vs 高意识型
es <- function(x, g1, g2) {
  a <- x[g1]; b <- x[g2]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  s <- sqrt(((length(a) - 1) * var(a) + (length(b) - 1) * var(b)) / (length(a) + length(b) - 2))
  (mean(a) - mean(b)) / s
}
gg <- d$alpha_group
cat(sprintf("\n  耗竭型 vs 高意识型 效应量（Cohen's d）：PHQ-9=%.2f｜GAD-7=%.2f｜疾病数=%.2f\n",
            es(d$phq9_total, gg == "耗竭型", gg == "高意识型"),
            es(d$gad7_total, gg == "耗竭型", gg == "高意识型"),
            es(d$disease_count, gg == "耗竭型", gg == "高意识型")))

if (length(krows)) save_tab(do.call(rbind, krows), "T23_known_groups")
cat("\n  ★ 判读要点：\n")
cat("    · PPV 强烈依赖患病率 —— 这正是「切点必须按场景标定」的实证依据；\n")
cat("      在同一 Sens/Spec 下，患病率 30% 与 5% 的 PPV 可相差数倍。\n")
cat("    · 已知组描述：分型与 PHQ-9 分级/疾病数的关联须结合分组定义审查；同源分组不构成独立效度证据。\n")
