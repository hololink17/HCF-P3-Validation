# =============================================================
# 10_analysis_18_anchor_vs_algorithm.R —— 计算型评分 vs 单题自评锚定（H2）
# -------------------------------------------------------------
# 对应蓝图 R2 的 H2：算法读数是否优于「每题一道自评锚定题」的朴素做法
# 仅 Prolific 完整版（四项金标准同时在测）
#
# 配对（金标准按理论方向统一为「越高越健康」）：
#   α₁ ↔ MAIA-2 三维度均分   （正）
#   α₂ ↔ AAQ-II              （负 → 取 −AAQ-II）
#   α₃ ↔ DCS                 （负 → 取 −DCS）
#   α₄ ↔ BADS 总分           （正）
#
# 每个 α 比较两个相关：
#   r_alg = r(α_i, 金标准_i)      算法读数
#   r_anc = r(锚定题_i, 金标准_i)  单题自评
#   r_xx  = r(α_i, 锚定题_i)       用于依赖相关检验
#
# 统计：Williams' t（1959；Steiger 1980）—— 两个「共享同一变量」的依赖相关之差
#   t = (r12 − r13)·sqrt{ (n−1)(1+r23) / [ 2·((n−1)/(n−3))·|R| + r̄²(1−r23)³ ] }
#   |R| = 1 − r12² − r13² − r23² + 2·r12·r13·r23 ；r̄ = (r12+r13)/2 ；df = n−3
# 辅助：bootstrap（5000 次）给出 Δr 的 95% CI
# 输出：T15_anchor_vs_algorithm.csv
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

# --- 本模块需要金标准量表数据（MAIA-2/AAQ-II/DCS/BADS），仅 Prolific 完整版队列施测 ---
# --- 故在见数运行模式下自动跳过，不影响见数侧任何结果 ---
if (identical(P3_COHORT, "jianshu")) {
  cat("\n（跳过：本模块仅适用于 Prolific 完整版队列——见数未施测金标准量表）\n")
} else {
GOLD_F <- base_path("p3_base_gold")
num <- c(ALPHA_AC, paste0("anchor", 1:4),
         "maia_three_dim", "aaq_total", "dcs_total", "bads_total")

g <- read_csv(GOLD_F, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
  mutate(across(any_of(num), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(maia_three_dim), !is.na(aaq_total), !is.na(dcs_total), !is.na(bads_total))
n_raw <- nrow(g)
if (DEDUP) {
  g <- g %>% arrange(participant_token, desc(data_source == "full")) %>%
    distinct(participant_token, .keep_all = TRUE)
}
cat(sprintf("\n=== 算法 vs 单题锚定（H2，完整版 %d → 去重后 %d）===\n", n_raw, nrow(g)))

# 公开 L2 不发布四个锚点题；此分析只能使用受控 L3，显式跳过而非报错。
if (!all(paste0("anchor", 1:4) %in% names(g))) {
  cat("  跳过：公开 L2 不含 anchor1–anchor4；本模块需要受控 L3。\n")
  .P3_MODULE_STATUS <- "skipped_controlled_l3_required"
  .P3_MODULE_REASON <- "需要受控 L3 的 anchor1–anchor4 锚点题"
} else {

# 金标准方向统一为「越高越健康」
g$gold1 <- g$maia_three_dim
g$gold2 <- -g$aaq_total
g$gold3 <- -g$dcs_total
g$gold4 <- g$bads_total
for (i in 1:4) g[[paste0("anc", i)]] <- g[[paste0("anchor", i)]]

williams_t <- function(r12, r13, r23, n) {
  det  <- 1 - r12^2 - r13^2 - r23^2 + 2 * r12 * r13 * r23
  rbar <- (r12 + r13) / 2
  denom <- 2 * ((n - 1) / (n - 3)) * det + rbar^2 * (1 - r23)^3
  if (!is.finite(denom) || denom <= 0) return(c(t = NA_real_, df = n - 3, p = NA_real_))
  t <- (r12 - r13) * sqrt((n - 1) * (1 + r23) / denom)
  c(t = t, df = n - 3, p = 2 * pt(-abs(t), n - 3))
}

GOLD_LAB <- c("MAIA-2 三维", "AAQ-II（已反向）", "DCS（已反向）", "BADS")
rows <- list()
set.seed(20260929)

for (i in 1:4) {
  a <- paste0("alpha", i); s <- paste0("anc", i); y <- paste0("gold", i)
  d <- g[, c(a, s, y)]; d <- d[complete.cases(d), ]
  n <- nrow(d)
  if (n < 20) next
  r_alg <- cor(d[[a]], d[[y]]); r_anc <- cor(d[[s]], d[[y]]); r_aa <- cor(d[[a]], d[[s]])
  w <- williams_t(r_alg, r_anc, r_aa, n)
  # bootstrap Δr
  B <- 5000; dr <- numeric(B)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    dr[b] <- cor(d[[a]][ii], d[[y]][ii]) - cor(d[[s]][ii], d[[y]][ii])
  }
  dr <- dr[is.finite(dr)]
  rows[[length(rows) + 1]] <- data.frame(
    alpha = a, gold = GOLD_LAB[i], n = n,
    r_algorithm = r_alg, r_anchor = r_anc, r_alg_anchor = r_aa,
    abs_r_algorithm = abs(r_alg), abs_r_anchor = abs(r_anc),
    delta_r = r_alg - r_anc,
    boot_lo = unname(quantile(dr, .025)), boot_hi = unname(quantile(dr, .975)),
    williams_t = w["t"], df = w["df"], p = w["p"],
    algorithm_better = ifelse(w["p"] < .05 & r_alg > r_anc, "是（显著）",
                       ifelse(r_alg > r_anc, "方向一致但未达显著", "否")),
    row.names = NULL)
  cat(sprintf("  %s ↔ %-14s n=%3d  r_算法=%+.3f  r_锚定=%+.3f  Δr=%+.3f [%+.3f,%+.3f]  t=%+.2f  p=%.4f  %s\n",
              a, GOLD_LAB[i], n, r_alg, r_anc, r_alg - r_anc,
              unname(quantile(dr, .025)), unname(quantile(dr, .975)),
              w["t"], w["p"], rows[[length(rows)]]$algorithm_better))
}

out <- do.call(rbind, rows)
save_tab(out, "T15_anchor_vs_algorithm")

cat("\n  ── 汇总 ──\n")
cat(sprintf("  4/4 方向一致（|r_算法| > |r_锚定|）：%s\n",
            ifelse(all(abs(out$r_algorithm) > abs(out$r_anchor)), "是", "否")))
cat(sprintf("  达显著（Williams t, p<.05）：%d / %d\n",
            sum(out$p < .05, na.rm = TRUE), nrow(out)))

cat("\n  ⚠️ 报告纪律：\n")
cat("     · 金标准已按理论方向统一（AAQ-II 与 DCS 取负），使「越高越健康」一致\n")
cat("     · 检验为依赖相关之差（Williams' t），非两独立相关比较\n")
cat("     · 锚定题与算法读数的相关（r_alg_anchor）一并报告，供读者判断二者是否同源\n")
cat("     · 若某行 r_锚定 与 r_算法 异号，应表述为「单题自评在方向上即不可靠」，\n")
cat("       而非仅「强度更弱」——这本身是结构读出优于单题自评的证据\n")

}  # 锚点列齐全时执行
}  # P3_COHORT != jianshu
