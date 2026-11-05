# =============================================================
# 10_analysis_21_acceptability.R —— 【零成本补强】可接受性 / 可行性 / 时间成本
# -------------------------------------------------------------
# 为什么要做：工具类论文（尤其 JMIR 系列）把「可用性与可行性」列为必备报告项；
#   此前 duration_seconds 与作答质量标记**从未被分析**（既有 19 个模块未覆盖）。
# 本模块回答四个问题：
#   ① 完成 20 题 HCF 需要多久？（中位/IQR，按队列与版本分层）
#   ② 换成 85 题金标准要多花多久？—— 这正是「20 题替代 85 题」的时间成本证据
#   ③ 作答质量如何？（无效标记率、缺失率、天花板/地板效应）
#   ④ 两队列的可行性是否一致？
# ⚠ 版本对比的性质：Prolific 60 题版（2026-05-31 投放）与 145 题版（2026-06-07~09 投放）
#   为**不同时点**的两批投放，属准实验对比，正文须声明这一混杂。
#   两版唯一差异 = 85 道金标准题（MAIA-2 37 + AAQ-II 7 + DCS 16 + BADS 25 = 85）。
# 口径：Prolific N=462（每人一行）｜见数 N=504
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

fmt_min <- function(s) ifelse(is.na(s), NA_character_, sprintf("%.1f", s / 60))

desc_dur <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 5) return(NULL)
  q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  data.frame(n = length(x), mean_s = round(mean(x), 1), sd_s = round(sd(x), 1),
             median_s = round(q[2], 1), q1_s = round(q[1], 1), q3_s = round(q[3], 1),
             iqr_s = round(q[3] - q[1], 1),
             min_s = round(min(x), 1), max_s = round(max(x), 1),
             median_min = round(q[2] / 60, 2), iqr_min = round((q[3] - q[1]) / 60, 2))
}

cat("\n=== 可接受性 / 可行性 ===\n")

# ---------- 一、当前队列的完成时长（按投放版本分层）----------
L1 <- load_cohorts()
d  <- L1[[1]]
coh <- d$cohort[1]

acc_rows <- list()
if ("duration_sec" %in% names(d) && any(is.finite(d$duration_sec))) {
  # 分版本（Prolific 有 full/short；见数只有一种源）
  srcs <- unique(na.omit(d$data_source))
  for (s in c(as.character(srcs), "ALL")) {
    sub <- if (s == "ALL") d else d[!is.na(d$data_source) & d$data_source == s, ]
    dd <- desc_dur(sub$duration_sec)
    if (is.null(dd)) next
    acc_rows[[length(acc_rows) + 1]] <- cbind(
      data.frame(cohort = coh, version = s),
      data.frame(n_items = ifelse(s == "full", 145,
                           ifelse(s == "short", 60, NA_real_)), stringsAsFactors = FALSE),
      dd)
    cat(sprintf("  %s｜%-6s 中位 %.1f 分钟（IQR %.1f–%.1f）｜均值 %.1f 分钟｜n=%d\n",
                coh, s, dd$median_min, dd$q1_s / 60, dd$q3_s / 60, dd$mean_s / 60, dd$n))
  }
} else {
  cat("  ⚠ 本队列无可用时长数据\n")
}
if (length(acc_rows)) save_tab(do.call(rbind, acc_rows), "T21_acceptability")

# ---------- 二、时间成本对比（跨队列，两侧写同一份）----------
LB <- load_both()
p  <- LB$prolific
j  <- LB$jianshu

get_med <- function(x) {
  x <- x[is.finite(x) & x > 0]
  if (length(x) < 5) return(NA_real_)
  as.numeric(median(x))
}
med_full  <- get_med(p$duration_sec[!is.na(p$data_source) & p$data_source == "full"])
med_short <- get_med(p$duration_sec[!is.na(p$data_source) & p$data_source == "short"])
med_js    <- get_med(j$duration_sec)
n_full    <- sum(!is.na(p$data_source) & p$data_source == "full" & is.finite(p$duration_sec))
n_short   <- sum(!is.na(p$data_source) & p$data_source == "short" & is.finite(p$duration_sec))

cost_rows <- list()
add_cost <- function(config, items, med_s, n, note) {
  cost_rows[[length(cost_rows) + 1]] <<- data.frame(
    config = config, n_items = items, n = n,
    median_s = round(med_s, 1), median_min = round(med_s / 60, 2),
    per_item_s = ifelse(is.na(items) | items == 0, NA_real_, round(med_s / items, 2)),
    note = note, stringsAsFactors = FALSE)
}
if (!is.na(med_short))
  add_cost("Prolific 60 题版（含 HCF 20 题）", 60, med_short, n_short, "实测")
if (!is.na(med_full))
  add_cost("Prolific 145 题版（含 HCF 20 题 + 金标准 85 题）", 145, med_full, n_full, "实测")
if (!is.na(med_js))
  add_cost("见数 60 题版（含 HCF 20 题）", 60, med_js, sum(is.finite(j$duration_sec)), "实测")

# 净增量：145 题版 − 60 题版 = 85 道金标准题的净时间成本
if (!is.na(med_full) && !is.na(med_short)) {
  delta <- med_full - med_short
  add_cost("★ 85 道金标准题的净时间成本（145 题版 − 60 题版）", 85, delta,
           min(n_full, n_short), "两版差值")
  # 按"每题均等"假设推算 20 题 HCF 的估计用时（须在正文声明为估算）
  hcf_est <- med_short * (20 / 60)
  add_cost("20 题 HCF 估计用时（按 60 题版均摊）", 20, hcf_est, n_short, "估算·假设每题均等")
  cat(sprintf("\n  ★ 时间成本核心结论：\n"))
  cat(sprintf("    60 题版（含 HCF 20 题）中位 %.1f 分钟 → 145 题版中位 %.1f 分钟\n",
              med_short / 60, med_full / 60))
  cat(sprintf("    两版中位时长的描述性差值为 %.1f 分钟；非同期版本对比存在队列与时间混杂\n",
              delta / 60))
  cat("    该时长差值不能证明筛查等效，也不能据此主张用一个量表替代另一个量表\n")
  # Mann–Whitney 检验（两版时长分布差异）
  a <- p$duration_sec[!is.na(p$data_source) & p$data_source == "full"]
  b <- p$duration_sec[!is.na(p$data_source) & p$data_source == "short"]
  wt <- suppressWarnings(wilcox.test(a, b))
  cat(sprintf("    两版时长 Mann–Whitney：W=%.0f，p=%.3e\n", wt$statistic, wt$p.value))
}

# 见数 vs Prolific 60 题版（同题量、跨队列的可行性对照）
if (!is.na(med_js) && !is.na(med_short)) {
  cat(sprintf("\n  跨队列可行性对照（同为 60 题版）：Prolific 中位 %.1f 分钟 vs 见数中位 %.1f 分钟（差 %.1f 分钟）\n",
              med_short / 60, med_js / 60, (med_js - med_short) / 60))
  add_cost("跨队列对照（60 题版）", 60, med_js, sum(is.finite(j$duration_sec)), "见数侧")
}
if (length(cost_rows)) save_tab(do.call(rbind, cost_rows), "T21_time_cost")

# ---------- 三、作答质量 ----------
q_rows <- list()
add_q <- function(cohort, item, value, note = "") {
  q_rows[[length(q_rows) + 1]] <<- data.frame(
    cohort = cohort, metric = item, value = value, note = note, stringsAsFactors = FALSE)
}
# 3.1 无效标记（Prolific 专有）
for (nm in c("prolific", "jianshu")) {
  dd <- LB[[nm]]
  if (all(is.na(dd$flag_attention)) && all(is.na(dd$flag_duration))) {
    add_q(dd$cohort[1], "作答质量校验字段", NA_real_, "该平台无此字段")
    next
  }
  add_q(dd$cohort[1], "注意力校验标记率",
        round(mean(dd$flag_attention == 1, na.rm = TRUE) * 100, 2), "% 被标记为无效")
  add_q(dd$cohort[1], "时长校验标记率",
        round(mean(dd$flag_duration == 1, na.rm = TRUE) * 100, 2), "% 被标记为无效")
}
# 3.2 关键变量缺失率
KEY <- c("alpha1", "alpha2", "alpha3", "alpha4", "phq9_total", "gad7_total", "ibs_positive")
for (nm in c("prolific", "jianshu")) {
  dd <- LB[[nm]]
  for (k in KEY) {
    if (!k %in% names(dd)) next
    add_q(dd$cohort[1], paste0("缺失率：", k),
          round(mean(is.na(dd[[k]])) * 100, 2), "%")
  }
}
# 3.3 天花板/地板（锚定题 6 级量表：最高档/最低档占比）
for (nm in c("prolific", "jianshu")) {
  dd <- LB[[nm]]
  for (a in c("anchor1", "anchor2", "anchor3", "anchor4")) {
    if (!a %in% names(dd)) next
    v <- suppressWarnings(as.numeric(dd[[a]]))
    v <- v[is.finite(v)]
    if (length(v) < 10) next
    add_q(dd$cohort[1], paste0("地板效应：", a, " = 最低档"),
          round(mean(v == min(v, na.rm = TRUE)) * 100, 2), "%")
    add_q(dd$cohort[1], paste0("天花板：", a, " = 最高档"),
          round(mean(v == max(v, na.rm = TRUE)) * 100, 2), "%")
  }
}
# 3.4 α 分数分布形态（偏度/峰度）
for (nm in c("prolific", "jianshu")) {
  dd <- LB[[nm]]
  for (a in c("alpha1", "alpha2", "alpha3", "alpha4")) {
    v <- suppressWarnings(as.numeric(dd[[a]])); v <- v[is.finite(v)]
    if (length(v) < 30) next
    s <- (sum((v - mean(v))^3) / length(v)) / (sd(v)^3)
    k <- (sum((v - mean(v))^4) / length(v)) / (sd(v)^4) - 3
    add_q(dd$cohort[1], paste0("偏度：", a), round(s, 3), "|偏度|<1 视为近似对称")
    add_q(dd$cohort[1], paste0("峰度（超额）：", a), round(k, 3), "接近 0 为近似正态")
  }
}
if (length(q_rows)) save_tab(do.call(rbind, q_rows), "T21_data_quality")
written21 <- c(if (length(acc_rows)) "T21_acceptability",
               if (length(cost_rows)) "T21_time_cost",
               if (length(q_rows)) "T21_data_quality")
cat("  → 实际写出表：", if (length(written21)) paste(written21, collapse = "、") else "无", "\n")
