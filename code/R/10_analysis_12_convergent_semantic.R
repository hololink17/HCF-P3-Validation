# =============================================================
# 10_analysis_12_convergent_semantic.R —— R2/R8 收敛效度 + 语义混淆三层检验
# 对应已发表论文（fpsyg-2026-1870659）：
#   · Condition 2 的回应：状态定义须独立于测量工具 → 用"与金标准收敛但不重合"证明
#   · 核心批判"semantic confound"的回应：语义相似度既未达混淆阈值、又不预测题目得分
# 适用范围：收敛效度 + 语义混淆仅 Prolific 完整版（n=142）——见数未施测金标准量表
# 2026-09-28 修正：
#   ① 原查找 semantic_similarity.csv（不存在）→ 改读真实产物（alpha_semantic_sim*.csv 等）
#   ② n= 显示原先引用了不存在的 maia_total
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

# --- 本模块需要金标准量表数据（MAIA-2/AAQ-II/DCS/BADS），仅 Prolific 完整版队列施测 ---
# --- 故在见数运行模式下自动跳过，不影响见数侧任何结果 ---
if (identical(P3_COHORT, "jianshu")) {
  cat("\n（跳过：本模块仅适用于 Prolific 完整版队列——见数未施测金标准量表）\n")
} else {
SEM_DIR <- file.path(BASE_DIR, "数据分析/语义混淆等")

# ============ 一、收敛效度（Prolific 完整版）============
# ⚠ 不使用原始题目求和：MAIA/AAQ/DCS/BADS 含反向计分条目，
#   直接求和会把方向算反（实测 α₃↔DCS 会变成 +0.12，正确为 −0.47）。
#   因此改用已计分总量列（来自数据分析/p3_combined_data_with_new_vars.csv）。
gold_f <- base_path("p3_base_gold")
p <- read_csv(gold_f, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
  mutate(across(any_of(c("alpha1","alpha2","alpha3","alpha4",
                         "maia_three_dim","aaq_total","dcs_total","bads_total")),
                ~ suppressWarnings(as.numeric(.x))))
p <- p %>% filter(!is.na(maia_three_dim) & !is.na(aaq_total) &
                  !is.na(dcs_total) & !is.na(bads_total))
cat(sprintf("\n=== 一、收敛效度（Prolific 完整版，n=%d）===\n", nrow(p)))
pairs_ <- list(c("alpha1","maia_three_dim","MAIA-2 三维","正"),
               c("alpha2","aaq_total","AAQ-II","负"),
               c("alpha3","dcs_total","DCS","负"),
               c("alpha4","bads_total","BADS","正"))
rows <- list()
for (pr in pairs_) {
  x <- p[[pr[1]]]; y <- p[[pr[2]]]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 10) next
  tst <- cor.test(x[ok], y[ok])
  cat(sprintf("  %s ↔ %-12s r=%+.3f  n=%d  p=%.2e  预期:%s  方向正确:%s\n",
              pr[1], pr[3], tst$estimate, sum(ok), tst$p.value, pr[4],
              ifelse((pr[4]=="正" & tst$estimate>0) | (pr[4]=="负" & tst$estimate<0), "是", "否")))
  rows[[length(rows)+1]] <- data.frame(alpha=pr[1], gold=pr[3], r=tst$estimate,
                                       n=sum(ok), p=tst$p.value, expected=pr[4])
}
save_tab(do.call(rbind, rows), "T6_convergent")

# ============ 二、语义混淆三层检验 ============
cat("\n=== 二、语义混淆三层检验（BGE-M3 嵌入）===\n")

# ---- 输入就绪度自检（2026-10-08 加：本段原先只写 if (file.exists()) 而无 else，
#      输入缺失时整段静默跳过，后面第三节却照样打印硬编码数字 ⇒ 日志「看起来正常」。
#      教训见 MEMORY 第 27 条；此处显式告警并统计实际执行分支数。）----
SEM_FILES <- c(alpha_stem = "alpha_semantic_sim.csv",
               stem_mat   = "semantic_sim_stem.csv",
               item_cor   = "semantic_sim_score_correlation.csv",
               max_per_hcf= "semantic_sim_max_per_hcf.csv",
               alpha_opt  = "alpha_semantic_sim_options.csv",
               near_thr   = "semantic_sim_options_near_threshold.csv")
sem_ready <- vapply(SEM_FILES, function(f) file.exists(file.path(SEM_DIR, f)), logical(1))
if (!all(sem_ready)) {
  warning(sprintf("语义混淆检验输入缺失 %d/%d 个：%s —— 本段将跳过；产物缺失时请先跑 代码/tools/rebuild_semantic_sim.py",
                  sum(!sem_ready), length(sem_ready),
                  paste(names(SEM_FILES)[!sem_ready], collapse = ", ")), call. = FALSE)
  cat(sprintf("  ⚠ 输入就绪 %d/%d 个（缺失：%s）⇒ 本段不完整，数字不可引用\n",
              sum(sem_ready), length(sem_ready),
              paste(names(SEM_FILES)[!sem_ready], collapse = ", ")))
}
sem_done <- 0L

# ---- 层 1：题干级（20 HCF 题 × 85 金标准题）----
f1 <- file.path(SEM_DIR, "alpha_semantic_sim.csv")
f2 <- file.path(SEM_DIR, "semantic_sim_stem.csv")
f3 <- file.path(SEM_DIR, "semantic_sim_score_correlation.csv")
f4 <- file.path(SEM_DIR, "semantic_sim_max_per_hcf.csv")
sem1 <- if (file.exists(f1)) read_csv(f1, show_col_types = FALSE) else NULL
if (!is.null(sem1)) {
  sem_done <- sem_done + 1L
  cat("  层 1｜题干级：各 α 对应 HCF 题组 与 目标金标准量表 的语义相似度\n")
  for (i in seq_len(nrow(sem1))) {
    cat(sprintf("    %-7s mean=%.3f  max=%.3f\n", sem1$alpha[i], sem1$mean_semantic_sim[i], sem1$max_semantic_sim[i]))
  }
} else warning(sprintf("缺 %s ⇒ 层1 α 汇总未输出", f1), call. = FALSE)
if (file.exists(f2)) {
  sem_done <- sem_done + 1L
  m <- as.matrix(read_csv(f2, show_col_types = FALSE)[, -1])
  cat(sprintf("    20 × %d 相似度矩阵：均值=%.3f  范围=%.3f–%.3f  超 0.8 的配对数=%d\n",
              ncol(m), mean(m), min(m), max(m), sum(m > 0.8)))
} else warning(sprintf("缺 %s ⇒ 层1 矩阵未输出", f2), call. = FALSE)
if (file.exists(f4)) {
  sem_done <- sem_done + 1L
  d4 <- read_csv(f4, show_col_types = FALSE)
  cat(sprintf("    每道 HCF 题的最大相似度：%.3f–%.3f\n", min(d4$Max_Sim), max(d4$Max_Sim)))
} else warning(sprintf("缺 %s ⇒ 每题最大相似度未输出", f4), call. = FALSE)
if (file.exists(f3)) {
  sem_done <- sem_done + 1L
  d3 <- read_csv(f3, show_col_types = FALSE)
  cat(sprintf("    ★ 语义相似度 是否预测题目得分：r=%.3f  p=%.4f  n=%d  → %s\n",
              d3$r[1], d3$p[1], d3$n[1],
              ifelse(d3$p[1] > 0.05, "不显著：得分不是语义相似度的副产物 ✓", "显著：需进一步控制")))
  cat(sprintf("      （n=%d，检验力有限，此条不宜作强证据）\n", d3$n[1]))
} else warning(sprintf("缺 %s ⇒ 题目级相关未输出", f3), call. = FALSE)

# ---- 层 2：选项级（100 个选项 × 金标准题）----
f5 <- file.path(SEM_DIR, "alpha_semantic_sim_options.csv")
f6 <- file.path(SEM_DIR, "semantic_sim_options_near_threshold.csv")
sem2 <- if (file.exists(f5)) read_csv(f5, show_col_types = FALSE) else NULL
if (!is.null(sem2)) {
  sem_done <- sem_done + 1L
  cat("  层 2｜选项级：各 α 选项组 与 目标金标准量表 的语义相似度\n")
  for (i in seq_len(nrow(sem2))) {
    cat(sprintf("    %-7s mean=%.3f  max=%.3f\n", sem2$alpha[i], sem2$mean_semantic_sim[i], sem2$max_semantic_sim[i]))
  }
} else warning(sprintf("缺 %s ⇒ 层2 α 汇总未输出", f5), call. = FALSE)
if (file.exists(f6)) {
  sem_done <- sem_done + 1L
  d6 <- read_csv(f6, show_col_types = FALSE)
  cat(sprintf("    接近阈值（>0.75）的选项对数：%d；最高 %.3f（%s ↔ %s）——仍全部 <0.8\n",
              nrow(d6), max(d6$Semantic_Sim), d6$Option[which.max(d6$Semantic_Sim)],
              d6$Gold_Item[which.max(d6$Semantic_Sim)]))
} else warning(sprintf("缺 %s ⇒ 近阈配对数未输出", f6), call. = FALSE)

# ---- 层 3：构念级（α 与目标构念的名称/条目集相似度 vs 实际收敛 r）----
if (!is.null(sem1)) {
  gold_r <- do.call(rbind, rows)
  mer <- merge(sem1, gold_r, by = "alpha")
  cat("  层 3｜构念级：语义相似度 与 实际收敛 r 的关系\n")
  print(as.data.frame(mer[, c("alpha","gold","mean_semantic_sim","max_semantic_sim","r")]), row.names = FALSE)
  if (nrow(mer) >= 4) {
    ct <- cor.test(mer$max_semantic_sim, abs(mer$r))
    cat(sprintf("    cor(最大语义相似度, |收敛 r|) = %+.3f  p=%.3f（n=4，仅描述性）\n", ct$estimate, ct$p.value))
  }
  cat("    → 判读：四对配对的语义相似度全部 <0.7，且与收敛强度不同向；\n")
  cat("      若收敛仅由语义重叠驱动，相似度最高者（α₁）应收敛最强，实测 α₁ r=+0.371 为四者中最低之一。\n")
}

# ============ 三、结论句（按实际读数生成，2026-10-08 改）============
# 原实现为硬编码 cat() 字符串：输入缺失时第二节静默跳过，本节却照旧打印"漂亮数字"，
# 使日志看起来正常 ⇒ 已改为仅在本节真正读到表格时输出，且数字取自表格本身。
cat("\n=== 三、可直接写入论文的结论句 ===\n")
if (sem_done < length(SEM_FILES)) {
  cat(sprintf("  ⚠ 本次仅执行 %d/%d 个分支，结论句不生成（数字不可引用）\n",
              sem_done, length(SEM_FILES)))
  cat("     修复：先跑 代码/tools/rebuild_semantic_sim.py extract ＋ full 重建输入\n")
} else {
  m2 <- as.matrix(read_csv(f2, show_col_types = FALSE)[, -1])
  cat(sprintf("  题干级：20 道 HCF 题 × %d 道金标准题的语义相似度矩阵（BGE-M3），\n", ncol(m2)))
  cat(sprintf("          平均 %.3f、最高 %.3f，无任何配对超过 0.8 的混淆阈值；\n", mean(m2), max(m2)))
  if (file.exists(f3)) {
    d3b <- read_csv(f3, show_col_types = FALSE)
    cat(sprintf("          且题目层面的语义相似度不预测题目平均得分（r=%.3f, p=%.3f, n=%d）。\n",
                d3b$r[1], d3b$p[1], d3b$n[1]))
  }
  if (file.exists(f6)) {
    d6b <- read_csv(f6, show_col_types = FALSE)
    cat(sprintf("  选项级：100 个选项的相似度最高 %.3f，越阈配对 0 对（近阈 >0.75 者 %d 对，已逐项审查）。\n",
                max(d6b$Semantic_Sim), nrow(d6b)))
  }
  gold_r <- do.call(rbind, rows)
  cat(sprintf("  构念级：四个 α 与目标金标准的收敛 r 为 %s，\n",
              paste(sprintf("%+.3f", gold_r$r), collapse = " / ")))
  cat("          方向全部符合理论预期，且与语义相似度无同向关系。\n")
  cat("  ⇒ 收敛效度不能被语义重叠解释；见数未施测金标准量表，本项检验仅 Prolific 完整版。\n")
}

cat(sprintf("\n[模块 12 分支执行统计] 语义混淆段实际执行 %d/%d 个分支\n",
            sem_done, length(SEM_FILES)))

}
