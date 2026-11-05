# =============================================================
# 10_analysis_11_cross_cohort.R —— R4.5 跨队列可复现性（论文最强段落）
# 本脚本用两个文化独立队列评价同一固定后台输出的跨队列统计表现
# 口径：Prolific 每人一行 n=462；见数 n=504
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_both(); P <- L$prolific; J <- L$jianshu   # 跨队列对照：强制读入两队列
AC <- c("alpha1","alpha2","alpha3","alpha4")
rows <- list()
cat("\n=== 跨队列：α → IBS（唯一口径统一后的对照）===\n")
for (a in AC) {
  vp <- auc_ci(P$ibs_positive, P[[a]]); vj <- auc_ci(J$ibs_positive, J[[a]])
  cat(sprintf("  %-7s  Prolific %.3f [%.3f,%.3f] (n=%.0f)   见数 %.3f [%.3f,%.3f] (n=%.0f)\n",
              a, vp["auc"], vp["lo"], vp["hi"], vp["n"], vj["auc"], vj["lo"], vj["hi"], vj["n"]))
  rows[[length(rows)+1]] <- data.frame(alpha=a, cohort="Prolific", AUC=vp["auc"], lo=vp["lo"], hi=vp["hi"], n=vp["n"])
  rows[[length(rows)+1]] <- data.frame(alpha=a, cohort="见数",     AUC=vj["auc"], lo=vj["lo"], hi=vj["hi"], n=vj["n"])
}
cvp <- cv_auc(P, AC); cvj <- cv_auc(J, AC)
cat(sprintf("  四α联合 5折CV:  Prolific %.3f±%.3f   见数 %.3f±%.3f\n", cvp["mean"], cvp["sd"], cvj["mean"], cvj["sd"]))

cat("\n=== ΔAUC 与 95% CI（bootstrap 3000 次）===\n")
for (a in AC) {
  b <- auc_diff_boot(P$ibs_positive, P[[a]], J$ibs_positive, J[[a]], B = 1000)
  cat(sprintf("  %-7s  ΔAUC = %+.3f  95%%CI [%+.3f, %+.3f]  → %s\n", a, b["delta"], b["lo"], b["hi"],
              ifelse(b["lo"] < 0 & b["hi"] > 0, "无显著差异（跨文化稳定）", "存在显著差异")))
  rows[[length(rows)+1]] <- data.frame(alpha=a, cohort="ΔAUC(Prolific-见数)",
                                       AUC=b["delta"], lo=b["lo"], hi=b["hi"], n=NA)
}
save_tab(do.call(rbind, rows), "T5_cross_cohort")

cat("\n=== 共同可比结局的两队列差异（效应量）===\n")
for (oc in c("phq9_total","gad7_total","ibs_positive","bmi","sleep_hour")) {
  if (!oc %in% names(P) || !oc %in% names(J)) next
  a <- P[[oc]]; b <- J[[oc]]
  tt <- tryCatch(t.test(a, b, var.equal = FALSE), error = function(e) NULL)
  if (is.null(tt)) next
  d <- (mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE))/sqrt((var(a,na.rm=TRUE)+var(b,na.rm=TRUE))/2)
  cat(sprintf("  %-14s Prolific %.2f  vs  见数 %.2f   Δ=%.2f  d=%.2f  p=%.2e\n",
              oc, mean(a,na.rm=TRUE), mean(b,na.rm=TRUE),
              mean(a,na.rm=TRUE)-mean(b,na.rm=TRUE), d, tt$p.value))
}
# 睡眠障碍：见数版该题选项结构与 Prolific 版不同（无「以上都没有」选项，阳性率恒为 100%），
# 按预案跨队列分析不纳入该结局。此处仅打印说明，不输出对比数字。
cat("  睡眠障碍      —— 见数版该题选项结构与 Prolific 版不同，按预案不纳入跨队列对比\n")

cat("\n=== 跨队列：分型与潜在痛苦（论文卖点）===\n")
cmp_typing <- function(d, nm) {
  n <- nrow(d)
  ag <- table(factor(d$alpha_group, levels = c("耗竭型","均衡型","高意识型","其他")))
  norm <- d %>% filter(phq9_total < 10)
  cat(sprintf("  %-9s 耗竭型 %.1f%% ｜ 均衡型 %.1f%% ｜ 高意识型 %.1f%% ｜ 其他 %.1f%%\n",
              nm, 100*ag[1]/n, 100*ag[2]/n, 100*ag[3]/n, 100*ag[4]/n))
  cat(sprintf("            PHQ-9 正常 n=%d → 潜在痛苦 %d (%.1f%%) ｜ 真正健康 %d (%.1f%%)\n",
              nrow(norm), sum(norm$potential_depression), 100*mean(norm$potential_depression),
              sum(norm$healthy_control), 100*mean(norm$healthy_control)))
  data.frame(cohort = nm, n = n,
             耗竭型 = 100*ag[1]/n, 均衡型 = 100*ag[2]/n, 高意识型 = 100*ag[3]/n, 其他 = 100*ag[4]/n,
             PHQ9正常_n = nrow(norm),
             潜在痛苦_n = sum(norm$potential_depression),
             潜在痛苦_pct = 100*mean(norm$potential_depression),
             真正健康_n = sum(norm$healthy_control),
             真正健康_pct = 100*mean(norm$healthy_control))
}
save_tab(rbind(cmp_typing(P, "Prolific"), cmp_typing(J, "见数")), "T5_cross_cohort_typing")

cat("\n  潜在痛苦 vs 真正健康 的 GAD-7 差（各队列内部）:\n")
for (nm in names(L)) {
  d <- L[[nm]] %>% filter(phq9_total < 10)
  a <- d$gad7_total[d$potential_depression]; b <- d$gad7_total[d$healthy_control]
  if (length(a) > 2 && length(b) > 2) {
    tt <- t.test(a, b)
    cat(sprintf("    %-9s 潜在痛苦 %.2f (n=%d) vs 真正健康 %.2f (n=%d)  Δ=%.2f  d=%.2f  p=%.2e\n",
                nm, mean(a), length(a), mean(b), length(b), mean(a)-mean(b),
                (mean(a)-mean(b))/sqrt((var(a)+var(b))/2), tt$p.value))
  }
}
