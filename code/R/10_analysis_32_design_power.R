# =============================================================
# 10_analysis_32_design_power.R —— C32 【设计层】E3a／E3b 最小样本量论证
# -------------------------------------------------------------
# 目的：为两项**待做采集**给出可复算的最小样本量，回应审稿人必问的
#       「你这个等价性／复现结论有多少把握？」
#
# 🔴 为什么报「设计层 MDE／最小 n」而不报「事后功效」：
#   事后功效是 p 值的确定性函数，不携带额外信息（Hoenig & Heisey 2001,
#   The American Statistician 55:19-24），会被直接判为统计误用。
#   正确做法是设计层反推：「要在给定把握下证明两套题等价，需多少样本？」
#   —— 与模块 14 同一口径（14 号管现有样本的 MDE，本模块管未来采集的 n）。
#
# 两项设计的定义（三类测量设计不可互替）：
#   · E3a 同批双题 = 同一批人 × 两套题目 × 同一时点
#         只证「测量工具层可迁移／收敛」；真正作用＝检验两套题的测量层是否
#         等价 ＋ 结果不利时提供归因分流。⚠️ 非必需项（E3b 不以它为前提）。
#   · E3b 独立样本施测新题 = 不同批人 × 新题
#         解除「跨域主张」禁令的唯一路径，可独立执行。
#
# 口径来源（全部沿用既有判据，不新造）：
#   ① 等价判据沿用模块 24 的 Romano 可忽略界 |Cliff's δ| < 0.147
#   ② δ→Cohen's d：d = √2·Φ⁻¹((1+δ)/2)
#   ③ TOST 单侧 α=0.05、功效 80% ⇒ (z₀.₉₅ + z₀.₈₀)²
#   ④ 配对设计 n = (z_a+z_b)²·2(1−ρ)/d²；ρ＝两套题同一人读数的相关（敏感性参数）
#   ⑤ 独立设计（原队列已存在）n_new = A/(1 − A/n_old)，A = (z_a+z_b)²/d²
#   ⑥ 精度口径：δ 的 95%CI 上界 < 可忽略界（Hanley–McNeil 近似，δ=0 最保守）
#
# ⚠️ 披露纪律：本模块**不做假设检验**，只做设计层样本量计算；
#   输出全部为聚合量（样本量、SD、判据阈值），符合 PUBLIC_BOUNDARY.md §2。
# =============================================================
# --- Windows 中文路径守卫 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

AC   <- c("alpha1", "alpha2", "alpha3", "alpha4")
ALLV <- c(AC, "H_score")

delta_to_d <- function(dl) sqrt(2) * qnorm((1 + dl) / 2)
Z2   <- (qnorm(0.95) + qnorm(0.80))^2          # TOST：单侧 α=.05、功效 .80
DELTAS <- c(0.200, 0.147, 0.100)               # 等价界（宽松 → 严格）
RHOS   <- c(0.50, 0.60, 0.70, 0.80, 0.85)      # 两套题同源相关（敏感性参数）

cat("=== C32 设计层：E3a／E3b 最小样本量论证 ===\n")
cat(sprintf("  TOST 口径：单侧 α=0.05、功效 80%% ⇒ (z+z)² = %.4f\n", Z2))
cat(sprintf("  δ→d 换算：%s\n",
            paste(sprintf("|δ|=%.3f→d=%.3f", DELTAS,
                          vapply(DELTAS, delta_to_d, 1)), collapse = " ｜ ")))

# ---------- 一、实测分布锚点（SD 尺度决定 margin 的绝对值）----------
L <- load_both()
cat("\n--- 一、实测分布锚点 ---\n")
anch <- list()
for (nm in c("prolific", "jianshu")) {
  d0 <- L[[nm]]
  for (a in ALLV) {
    if (!a %in% names(d0)) next
    x <- suppressWarnings(as.numeric(d0[[a]])); x <- x[is.finite(x)]
    anch[[length(anch) + 1]] <- data.frame(
      cohort = d0$cohort[1], alpha = a, n = length(x),
      mean = round(mean(x), 4), sd = round(sd(x), 4),
      lo = round(min(x), 4), hi = round(max(x), 4), stringsAsFactors = FALSE)
    cat(sprintf("  %-9s %-8s n=%4d  mean=%8.4f  sd=%7.4f\n",
                d0$cohort[1], a, length(x), mean(x), sd(x)))
  }
}
anch_tab <- do.call(rbind, anch)
save_tab(anch_tab, "T32_design_anchor")
# 合并 SD（两队列）作为 σ 的锚
pool_sd <- sapply(ALLV, function(a) {
  v <- c()
  for (nm in c("prolific", "jianshu")) {
    if (!a %in% names(L[[nm]])) next
    x <- suppressWarnings(as.numeric(L[[nm]][[a]])); v <- c(v, x[is.finite(x)])
  }
  sd(v)
})
cat("  合并 SD：", paste(sprintf("%s=%.4f", names(pool_sd), pool_sd), collapse = "  "), "\n")

# ---------- 二、E3a：同批双题（配对 TOST）----------
cat("\n--- 二、E3a 同批双题：所需【配对完成人数】---\n")
cat("  n = (z_a+z_b)²·2(1−ρ)/d² ；ρ 越小（两套题越独立）需求越大\n")
rowsA <- list()
for (rho in RHOS) for (dl in DELTAS) {
  d <- delta_to_d(dl)
  n <- ceiling(Z2 * 2 * (1 - rho) / d^2)
  rowsA[[length(rowsA) + 1]] <- data.frame(
    rho = rho, delta_margin = dl, cohen_d = round(d, 4), n_pairs = n,
    stringsAsFactors = FALSE)
}
tabA <- do.call(rbind, rowsA)
cat(sprintf("  %-6s ｜ %-16s %-16s %-16s\n", "ρ",
            sprintf("|δ|<%.3f", DELTAS[1]), sprintf("|δ|<%.3f", DELTAS[2]),
            sprintf("|δ|<%.3f", DELTAS[3])))
for (rho in RHOS) {
  sub <- tabA[tabA$rho == rho, ]
  cat(sprintf("  %.2f   ｜ %-16d %-16d %-16d\n", rho, sub$n_pairs[1],
              sub$n_pairs[2], sub$n_pairs[3]))
}
save_tab(tabA, "T32_e3a_paired_n")

# ---------- 三、E3b：独立样本施测新题（不平衡两独立组 TOST）----------
cat("\n--- 三、E3b 独立样本：原队列已存在时所需【新招人数】---\n")
n_old_map <- c(Prolific = nrow(L$prolific), jianshu = nrow(L$jianshu))
rowsB <- list()
for (k in seq_along(n_old_map)) {
  n_old <- n_old_map[[k]]
  for (dl in DELTAS) {
    d <- delta_to_d(dl); A <- Z2 / d^2
    n_eq  <- ceiling(2 * A)
    n_new <- ceiling(A / (1 - A / n_old))
    rowsB[[length(rowsB) + 1]] <- data.frame(
      anchor_cohort = names(n_old_map)[k], n_old = n_old,
      delta_margin = dl, cohen_d = round(d, 4),
      n_per_group_balanced = n_eq, n_new_required = n_new,
      stringsAsFactors = FALSE)
    cat(sprintf("  已知 %-9s n=%3d ｜ |δ|<%.3f ｜ 等组需 %4d/组 ｜ 新招 %4d\n",
                names(n_old_map)[k], n_old, dl, n_eq, n_new))
  }
}
tabB <- do.call(rbind, rowsB)
save_tab(tabB, "T32_e3b_independent_n")

# ---------- 四、精度口径（δ 的 95%CI 上界 < 可忽略界）----------
cat("\n--- 四、精度口径：δ 的 95%CI 上界 < 可忽略界（Hanley–McNeil，δ=0 最保守）---\n")
se_delta0 <- function(n1, n2) 2 * sqrt((0.25 + 0.0833 * (n1 + n2 - 2)) / (n1 * n2))
rowsP <- list()
for (dl in DELTAS) {
  f <- function(n) 1.96 * se_delta0(n, n) - dl
  n_eq <- ceiling(uniroot(f, c(10, 200000))$root)
  for (k in seq_along(n_old_map)) {
    n_old <- n_old_map[[k]]
    g <- function(n) 1.96 * se_delta0(n, n_old) - dl
    n_new <- if (g(2) < 0) 2 else ceiling(uniroot(g, c(2, 500000))$root)
    rowsP[[length(rowsP) + 1]] <- data.frame(
      delta_margin = dl, n_per_group_balanced = n_eq,
      anchor_cohort = names(n_old_map)[k], n_old = n_old,
      n_new_required = n_new, stringsAsFactors = FALSE)
    cat(sprintf("  |δ|<%.3f：等组需 %5d/组 ｜ 已知 %s n=%d 时新招 %5d\n",
                dl, n_eq, names(n_old_map)[k], n_old, n_new))
  }
}
save_tab(do.call(rbind, rowsP), "T32_precision_delta_ci")

# ---------- 五、E3a 邀请量（含重测留存率余量）----------
cat("\n--- 五、E3a 邀请量（含重测留存率余量）---\n")
rowsR <- list()
for (tgt in c(87, 100, 120, 160)) for (ret in c(0.80, 0.60, 0.50)) {
  rowsR[[length(rowsR) + 1]] <- data.frame(
    target_pairs = tgt, retention = ret,
    invite_required = ceiling(tgt / ret), stringsAsFactors = FALSE)
  cat(sprintf("  目标 %3d 对 ｜ 留存 %.0f%% ⇒ 需邀请 %4d 人\n",
              tgt, ret * 100, ceiling(tgt / ret)))
}
save_tab(do.call(rbind, rowsR), "T32_e3a_invite_size")

# ---------- 六、结论口径 ----------
cat("\n--- 六、设计层结论（供论文 Methods 引用）---\n")
n_e3a_ref <- tabA$n_pairs[tabA$rho == 0.70 & tabA$delta_margin == 0.147]
n_e3b_ref <- tabB$n_new_required[tabB$delta_margin == 0.147 &
                                   tabB$anchor_cohort == "Prolific"]
cat(sprintf("  · E3a（ρ=0.70、|δ|<0.147）最小配对完成数：%d 对\n", n_e3a_ref))
cat(sprintf("  · E3b（|δ|<0.147、以 Prolific 为参照）最小新招数：%d 人\n", n_e3b_ref))
cat(sprintf("  · 稳健取值建议：E3a 目标 100 对（留存 60%% ⇒ 邀请 %d 人）；\n",
            ceiling(100 / 0.60)))
cat("     E3b 目标 300 人（覆盖结构复现所需量级，远超等价检验下限）。\n")
cat("  · ⚠️ 判读纪律：E3b 除「分布等价」外还须复现核心结构，\n")
cat("     故实际样本量应取「等价下限」与「结构复现所需量级」之大者。\n")
