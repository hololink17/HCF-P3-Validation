# =============================================================
# 10_analysis_24_item_equivalence.R —— R2 【零成本补强】题目等价性与跨队列分布等价性
# -------------------------------------------------------------
# 为什么这样做（而不是 CFA 测量不变性）：
#   HCF 四 α 是 **formative/computational index**（形成性指数），不是反身性潜变量，
#   传统测量不变性（configural / metric / scalar MI）对形成性指数在测量学上不适用。
#   本模块给出**形成性框架下的等价性替代证据**。
#
# 公开 L2 不提供内部题目作答，因此公开运行只执行 α/H 的跨队列分布比较；
# 需要受控题目级输入的分支必须显式跳过。
#
# 本模块给出两层等价性证据：
#   ① 【题目级】Prolific 内部跨版本（145 题完整版 vs 60 题精简版，同 20 题分值口径一致）
#      —— 检验题目行为的时间稳定性（两版为 2026-05-31 与 2026-06-07~09 两批投放）
#   ② 【α 层面】跨队列分布等价性（Prolific ↔ 见数）
#      —— 检验同一算法口径在两个独立样本上的分数分布是否一致
# 本模块为跨队列分析：无论 P3_COHORT 为哪一种，两侧都写出同一份结果（与 T5/T10/T11 同规）。
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

ITEM_COLS <- paste0("hcf_", 5:24)
ALPHA_COLS <- c("alpha1", "alpha2", "alpha3", "alpha4")

# Cliff's delta：δ = 2U/(n1·n2) − 1 ∈ [−1, 1]（0 = 无差异）
cliff_delta <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  n1 <- length(x); n2 <- length(y)
  if (n1 < 10 || n2 < 10) return(NA_real_)
  r <- rank(c(x, y))
  U <- sum(r[seq_len(n1)]) - n1 * (n1 + 1) / 2
  2 * U / (n1 * n2) - 1
}
# 效应量判读（Romano et al. 2006）：|δ|<0.147 可忽略；<0.33 小；<0.474 中；≥0.474 大
delta_label <- function(d) {
  ifelse(is.na(d), NA_character_,
    ifelse(abs(d) < 0.147, "可忽略", ifelse(abs(d) < 0.33, "小",
      ifelse(abs(d) < 0.474, "中", "大"))))
}
# 重叠系数 OVL（基于核密度估计，0–1，1 = 完全重叠）
ovl <- function(x, y) {
  x <- x[is.finite(x)]; y <- y[is.finite(y)]
  if (length(x) < 20 || length(y) < 20) return(NA_real_)
  rng <- range(c(x, y)); g <- seq(rng[1], rng[2], length.out = 512)
  dx <- density(x, from = rng[1], to = rng[2], n = 512)$y
  dy <- density(y, from = rng[1], to = rng[2], n = 512)$y
  sum(pmin(dx, dy)) * (g[2] - g[1])
}

cat("\n=== 题目等价性与跨队列分布等价性（形成性框架）===\n")

# ---------- 一、题目级：Prolific 跨版本（同 20 题，分值口径一致）----------
pk <- read_csv(base_path("p3_base_prolific"),
               col_types = cols(.default = col_character()), show_col_types = FALSE)
has_items <- all(ITEM_COLS %in% names(pk))
if (has_items) {
  pk <- pk %>% mutate(across(all_of(ITEM_COLS), ~ suppressWarnings(as.numeric(.x))))
} else {
  cat("  公开 L2 不含 hcf_5–hcf_24 原始题项；题目级部分需要受控 L3，仅继续执行 α 层面分析。\n")
}

# ⚠ 口径修正（2026-10-02）：prolific_clean.csv 属「清洗层」= 471 行（含 9 行重复），
#   而全文主口径为「每人一行」= 462。此前本模块直接读清洗层，导致
#   完整版 n=142 + 精简版 n=329 = 471，与主口径 462 不一致（差 9 行）。
#   此处按 10_analysis_00_common.R 的同一规则去重：同一 participant_token
#   既参加完整版又参加精简版时，优先保留 full。
#   → 修正后：完整版 142 + 精简版 320 = 462，与全文口径一致。
if ("participant_token" %in% names(pk) && "data_source" %in% names(pk)) {
  n_before <- nrow(pk)
  pk <- pk %>% arrange(participant_token, desc(data_source == "full")) %>%
    distinct(participant_token, .keep_all = TRUE)
  cat(sprintf("  [口径] prolific_clean.csv %d 行 → 去重后 %d 行（每人一行，重复者优先保留 full）\n",
              n_before, nrow(pk)))
}

has_src <- has_items && "data_source" %in% names(pk)
rows <- list()
if (has_src) {
  pf <- pk[!is.na(pk$data_source) & pk$data_source == "full", ]
  ps <- pk[!is.na(pk$data_source) & pk$data_source == "short", ]
  cat(sprintf("  题目级对比：Prolific 完整版 n=%d vs 精简版 n=%d（同 20 题分值口径）\n",
              nrow(pf), nrow(ps)))
  for (it in ITEM_COLS) {
    xf <- pf[[it]]; xs <- ps[[it]]
    if (sum(is.finite(xf)) < 30 || sum(is.finite(xs)) < 30) next
    ks <- suppressWarnings(ks.test(xf, xs))
    dl <- cliff_delta(xf, xs)
    rows[[length(rows) + 1]] <- data.frame(
      scope = "Prolific 完整版 vs 精简版（同 20 题）",
      item = it,
      n_a = sum(is.finite(xf)), n_b = sum(is.finite(xs)),
      mean_a = round(mean(xf, na.rm = TRUE), 3), mean_b = round(mean(xs, na.rm = TRUE), 3),
      sd_a = round(sd(xf, na.rm = TRUE), 3), sd_b = round(sd(xs, na.rm = TRUE), 3),
      ks_stat = round(as.numeric(ks$statistic), 4),
      ks_p = format.pval(ks$p.value, digits = 3),
      cliff_delta = round(dl, 4), effect_size = delta_label(dl),
      ovl = round(ovl(xf, xs), 3),
      stringsAsFactors = FALSE)
  }
} else {
  cat("  跳过题目级跨版本对比（需要受控 L3 原始题项及 data_source 版本标记）。\n")
}

# ---------- 二、α 层面：跨队列分布等价性 ----------
LB <- load_both()
p <- LB$prolific; j <- LB$jianshu
cat(sprintf("\n  α 层面跨队列对比：Prolific n=%d（%s） vs 见数 n=%d（%s）\n",
            nrow(p), p$cohort[1], nrow(j), j$cohort[1]))

# 2.1 各 α 与 H 的分布等价性
for (a in c(ALPHA_COLS, "H_score")) {
  xa <- suppressWarnings(as.numeric(p[[a]])); xb <- suppressWarnings(as.numeric(j[[a]]))
  if (sum(is.finite(xa)) < 30 || sum(is.finite(xb)) < 30) next
  ks <- suppressWarnings(ks.test(xa, xb))
  dl <- cliff_delta(xa, xb)
  ov <- ovl(xa, xb)
  rows[[length(rows) + 1]] <- data.frame(
    scope = "跨队列 α 层面（Prolific vs 见数）",
    item = a,
    n_a = sum(is.finite(xa)), n_b = sum(is.finite(xb)),
    mean_a = round(mean(xa, na.rm = TRUE), 3), mean_b = round(mean(xb, na.rm = TRUE), 3),
    sd_a = round(sd(xa, na.rm = TRUE), 3), sd_b = round(sd(xb, na.rm = TRUE), 3),
    ks_stat = round(as.numeric(ks$statistic), 4),
    ks_p = format.pval(ks$p.value, digits = 3),
    cliff_delta = round(dl, 4), effect_size = delta_label(dl),
    ovl = round(ov, 3),
    stringsAsFactors = FALSE)
  cat(sprintf("    %-8s 均值 %.3f vs %.3f｜Cliff's δ=%+.3f（%s）｜重叠系数=%.3f\n",
              a, mean(xa, na.rm = TRUE), mean(xb, na.rm = TRUE), dl, delta_label(dl), ov))
}

# 2.2 分型构成比例一致性（P33/P67 各自分档，故比较的是比例结构）
gp <- prop.table(table(p$alpha_group)); gj <- prop.table(table(j$alpha_group))
lvl <- union(names(gp), names(gj))
cmp <- data.frame(
  alpha_group = lvl,
  prop_prolific = round(as.numeric(gp[lvl]), 4),
  prop_jianshu  = round(as.numeric(gj[lvl]), 4),
  stringsAsFactors = FALSE)
cmp$diff <- round(cmp$prop_prolific - cmp$prop_jianshu, 4)
cat("\n  分型构成比例（各队列自行按 P33/P67 分档，比较的是比例结构）\n")
print(cmp, row.names = FALSE)
save_tab(cmp, "T24_alpha_typing_proportions")

if (length(rows)) {
  tab <- do.call(rbind, rows)
  save_tab(tab, "T24_item_equivalence")

  cat("\n  ★ 整体判据：\n")
  it_tab <- tab[grepl("完整版 vs 精简版", tab$scope), ]
  if (nrow(it_tab)) {
    cat(sprintf("    · 【题目级】Prolific 跨版本：|δ| 中位 %.3f｜最大 %.3f｜判定「可忽略」%d/%d 题｜重叠系数中位 %.3f\n",
                median(abs(it_tab$cliff_delta), na.rm = TRUE),
                max(abs(it_tab$cliff_delta), na.rm = TRUE),
                sum(abs(it_tab$cliff_delta) < 0.147, na.rm = TRUE), nrow(it_tab),
                median(it_tab$ovl, na.rm = TRUE)))
  }
  al_tab <- tab[grepl("跨队列", tab$scope), ]
  if (nrow(al_tab)) {
    cat(sprintf("    · 【α 层面】跨队列：|δ| 中位 %.3f｜最大 %.3f（%s）｜重叠系数中位 %.3f\n",
                median(abs(al_tab$cliff_delta), na.rm = TRUE),
                max(abs(al_tab$cliff_delta), na.rm = TRUE),
                al_tab$item[which.max(abs(al_tab$cliff_delta))],
                median(al_tab$ovl, na.rm = TRUE)))
  }
  cat("    · 判读：两个大样本下 KS 检验对微小差异也会显著，故应以**效应量（Cliff's δ）\n")
  cat("      与重叠系数**判读分布差异；|δ|<0.147 视为可忽略。\n")
  cat("    · 结论写法建议：形成性指数的等价性以「效应量可忽略 + 分布高度重叠」为准，\n")
  cat("      而非传统 MI 的 ΔCFI / ΔRMSEA（后者对形成性指数不适用）。\n")
  cat("    · 公开输入不包含内部题目级数据，故仅报告 α/H 层面的分布等价性。\n")
}
