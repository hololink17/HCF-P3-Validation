# =============================================================
# 10_analysis_01_descriptives.R —— R1 样本、分型与描述统计（TABLE 1）
# 对应已发表论文：Condition 2（状态定义独立于工具）——先证明"状态可被观察"
# 两队列同规则分型（profile 型，见 R/00_shared_rules.R::add_alpha_typing）
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts()
res <- list(); typ <- list(); flow <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  res[[nm]] <- data.frame(
    cohort = nm, N = nrow(d),
    age_18_35 = mean(d$age_group %in% "18-35", na.rm=TRUE)*100,
    age_36_45 = mean(d$age_group %in% "36-45", na.rm=TRUE)*100,
    age_46_60 = mean(d$age_group %in% "46-60", na.rm=TRUE)*100,
    age_60p   = mean(d$age_group %in% "60+",  na.rm=TRUE)*100,
    female_pct = mean(d$gender %in% "女", na.rm=TRUE)*100,
    phq9_mean = mean(d$phq9_total, na.rm=TRUE), phq9_sd = sd(d$phq9_total, na.rm=TRUE),
    gad7_mean = mean(d$gad7_total, na.rm=TRUE), gad7_sd = sd(d$gad7_total, na.rm=TRUE),
    depression_pct = mean(d$phq9_total >= 10, na.rm=TRUE)*100,
    ibs_pct = mean(d$ibs_positive, na.rm=TRUE)*100, ibs_n = sum(d$ibs_positive, na.rm=TRUE),
    bmi_mean = mean(d$bmi, na.rm=TRUE),
    a1 = mean(d$alpha1, na.rm=TRUE), a2 = mean(d$alpha2, na.rm=TRUE),
    a3 = mean(d$alpha3, na.rm=TRUE), a4 = mean(d$alpha4, na.rm=TRUE),
    H_mean = mean(d$H_score, na.rm=TRUE)
  )
  cat(sprintf("\n--- %s (N=%d) ---\n", nm, nrow(d)))
  cat(sprintf("  PHQ-9 %.1f±%.1f | GAD-7 %.1f±%.1f | 抑郁阳性 %.1f%% | IBS %.1f%% (%d)\n",
              res[[nm]]$phq9_mean, res[[nm]]$phq9_sd, res[[nm]]$gad7_mean, res[[nm]]$gad7_sd,
              res[[nm]]$depression_pct, res[[nm]]$ibs_pct, res[[nm]]$ibs_n))
  cat(sprintf("  BMI %.1f | α: %.3f / %.3f / %.3f / %.3f | H %.3f\n",
              res[[nm]]$bmi_mean, res[[nm]]$a1, res[[nm]]$a2, res[[nm]]$a3, res[[nm]]$a4, res[[nm]]$H_mean))
  cat(sprintf("  年龄 18-35 %.1f%% / 36-45 %.1f%% / 46-60 %.1f%% / 60+ %.1f%% ｜ 女性 %.1f%%\n",
              res[[nm]]$age_18_35, res[[nm]]$age_36_45, res[[nm]]$age_46_60, res[[nm]]$age_60p, res[[nm]]$female_pct))

  # ---- 分型分布（profile 型；含「其他」）----
  tb <- table(factor(d$alpha_group, levels = c("耗竭型","均衡型","高意识型","其他")))
  cat("  α 分型:")
  for (k in names(tb)) cat(sprintf("  %s %d (%.1f%%)", k, tb[k], 100*tb[k]/nrow(d)))
  cat("\n")
  norm <- d %>% filter(phq9_total < 10)
  cat(sprintf("  PHQ-9 正常者 n=%d ｜ 潜在痛苦（耗竭型）%d (%.1f%%) ｜ 真正健康（高意识型）%d (%.1f%%)\n",
              nrow(norm), sum(norm$potential_depression), 100*mean(norm$potential_depression),
              sum(norm$healthy_control), 100*mean(norm$healthy_control)))
  typ[[nm]] <- data.frame(cohort = nm, alpha_group = names(tb), n = as.integer(tb),
                          pct = as.numeric(100*tb/nrow(d)),
                          potential_depression_n = sum(norm$potential_depression),
                          phq9_normal_n = nrow(norm))

  # ---- 样本构成（2026-10-05 新增）----
  # 论文 CONSORT 式流程图与 Methods 的唯一引用值：
  #   Prolific N=462 = 142 完整版（145 题）+ 320 精简版（60 题）
  #   见数     N=504（后台两批）
  # ⚠️ 对外禁用 368 / 471 / 519 作为分析样本量（471 是清洗层未去重行数，519 是提交总数）。
  ds <- table(factor(d$data_source, levels = c("full", "short",
                                               "cleaned_60q.csv", "cleaned_60q_merged.csv")))
  ds <- ds[ds > 0]
  cat(sprintf("  样本构成：%s ｜ 合计 %d\n",
              paste(sprintf("%s %d", names(ds), as.integer(ds)), collapse = " + "), nrow(d)))
  # 🔴 见数队列为单一问卷版本（无完整版／精简版之分）——此时 n_full/n_short 必须写 NA
  #    而非 0，否则会被误读为「该队列没有完整版样本」（0 与「不适用」语义不同）。
  n_full <- sum(d$data_source == "full", na.rm = TRUE)
  n_short <- sum(d$data_source == "short", na.rm = TRUE)
  has_split <- (n_full + n_short) > 0
  comp <- if (has_split) {
    paste(sprintf("%s=%d", names(ds), as.integer(ds)), collapse = " | ")
  } else if (length(ds) > 0) {
    paste0("single-version cohort (no full/short split) | source=",
           paste(sprintf("%s=%d", names(ds), as.integer(ds)), collapse = " + "))
  } else {
    # 无 full/short 且 data_source 无匹配级别 —— 不留空的 "source=" 尾巴
    "single-version cohort (no full/short split)"
  }
  flow[[length(flow) + 1]] <- data.frame(
    cohort = nm, N_analysis = nrow(d),
    n_full = if (has_split) n_full else NA_integer_,
    n_short = if (has_split) n_short else NA_integer_,
    composition = comp,
    stringsAsFactors = FALSE)
}
save_tab(do.call(rbind, res), "T1_descriptives")
save_tab(do.call(rbind, typ), "T1_alpha_typing")
save_tab(do.call(rbind, flow), "T1_sample_flow")
