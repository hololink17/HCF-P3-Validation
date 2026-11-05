# =============================================================
# 10_analysis_22_incremental_mtmm.R —— R2/R8 【零成本补强】增量效度 + MTMM 收敛-区分矩阵
# -------------------------------------------------------------
# 为什么要做（回应最致命的一句质疑）：
#   审稿人最可能的攻击是「为什么不直接用现成的 MAIA-2 / AAQ-II / DCS / BADS？」
#   模块 12 的收敛效度（T6，r = .371 ~ .647）只证明"所测东西有关联"，
#   **不证明 α 提供了这四个成熟量表之外的新信息**。
#   本模块用两条证据回应：
#     ① 增量效度：控制四个金标准后，α 是否仍显著提升对 PHQ-9 / GAD-7 的解释与判别
#        （连续结局看 ΔR²、二分类结局看 ΔAUC + DeLong）
#     ② MTMM 式收敛-区分矩阵：补齐 4×4 全矩阵，
#        验证"目标配对的 |r| 应高于非目标配对"（否则只是共同方法变异）
# 适用范围：**仅 Prolific 完整版（n=142）**——见数未施测金标准量表，故见数模式自动跳过。
# 数据源：数据分析/p3_combined_data_with_new_vars.csv
#   ⚠ 不使用原始题目求和：金标准量表含反向计分条目，直接求和会把方向算反
#     （实测 α₃↔DCS 会变成 +0.12，正确为 −0.47）。故一律使用已计分总量列。
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

if (identical(P3_COHORT, "jianshu")) {
  cat("\n（跳过：本模块仅适用于 Prolific 完整版队列——见数未施测金标准量表）\n")
} else {

GOLD_F <- base_path("p3_base_gold")
stopifnot(file.exists(GOLD_F))
p <- read_csv(GOLD_F, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
  mutate(across(any_of(c("alpha1", "alpha2", "alpha3", "alpha4",
                         "maia_three_dim", "aaq_total", "dcs_total", "bads_total",
                         "phq9_total", "gad7_total")),
                ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(maia_three_dim) & !is.na(aaq_total) & !is.na(dcs_total) & !is.na(bads_total))
p$phq9_pos <- as.numeric(p$phq9_total >= 10)
p$gad7_pos <- as.numeric(p$gad7_total >= 10)

cat(sprintf("\n=== 增量效度 + MTMM（Prolific 完整版，n=%d）===\n", nrow(p)))

GOLDS <- c(maia_three_dim = "MAIA-2 三维", aaq_total = "AAQ-II",
           dcs_total = "DCS", bads_total = "BADS")
ALPHAS <- c(alpha1 = "α₁", alpha2 = "α₂", alpha3 = "α₃", alpha4 = "α₄")
# 理论预期的目标配对（其余为非目标配对）
TARGET <- c(alpha1 = "maia_three_dim", alpha2 = "aaq_total",
            alpha3 = "dcs_total",       alpha4 = "bads_total")

# ============ 一、MTMM 收敛-区分矩阵（4 × 4 = 16 对）============
cat("\n=== 一、MTMM 收敛-区分矩阵 ===\n")
mrows <- list()
for (a in names(ALPHAS)) {
  for (g in names(GOLDS)) {
    x <- p[[a]]; y <- p[[g]]
    ok <- is.finite(x) & is.finite(y)
    if (sum(ok) < 10) next
    tst <- suppressWarnings(cor.test(x[ok], y[ok]))
    is_t <- identical(unname(TARGET[a]), g)
    mrows[[length(mrows) + 1]] <- data.frame(
      alpha = ALPHAS[a], gold = GOLDS[g], r = round(as.numeric(tst$estimate), 4),
      n = sum(ok), p = format.pval(tst$p.value, digits = 3),
      pair_type = ifelse(is_t, "目标配对（收敛）", "非目标配对（区分）"),
      stringsAsFactors = FALSE)
  }
}
mtmm <- do.call(rbind, mrows)
save_tab(mtmm, "T22_mtmm")

# 收敛 vs 区分 的强度对照（判别效度的核心指标）
tgt <- mtmm$r[mtmm$pair_type == "目标配对（收敛）"]
non <- mtmm$r[mtmm$pair_type == "非目标配对（区分）"]
cat(sprintf("  目标配对 |r| 范围 %.3f–%.3f（均值 %.3f，n=%d）\n",
            min(abs(tgt)), max(abs(tgt)), mean(abs(tgt)), length(tgt)))
cat(sprintf("  非目标配对 |r| 范围 %.3f–%.3f（均值 %.3f，n=%d）\n",
            min(abs(non)), max(abs(non)), mean(abs(non)), length(non)))
cat(sprintf("  ★ 判别效度判据：目标配对强度是否系统性高于非目标配对 —— %s\n",
            ifelse(mean(abs(tgt)) > mean(abs(non)),
                   sprintf("是（%.3f > %.3f）", mean(abs(tgt)), mean(abs(non))),
                   sprintf("否（%.3f ≤ %.3f）", mean(abs(tgt)), mean(abs(non))))))

# ============ 二、增量效度 ============
cat("\n=== 二、增量效度（控制四个金标准后，α 是否仍有增量）===\n")
GOLD_VARS <- c("maia_three_dim", "aaq_total", "dcs_total", "bads_total")
ALPHA_VARS <- c("alpha1", "alpha2", "alpha3", "alpha4")
irows <- list()

add_inc <- function(outcome, otype, m0_r2, m1_r2, extra) {
  irows[[length(irows) + 1]] <<- cbind(
    data.frame(outcome = outcome, type = otype,
               R2_gold_only = round(m0_r2, 4), R2_gold_plus_alpha = round(m1_r2, 4),
               delta_R2 = round(m1_r2 - m0_r2, 4),
               delta_R2_pct = round((m1_r2 - m0_r2) / max(m0_r2, 1e-9) * 100, 1),
               stringsAsFactors = FALSE),
    extra)
}

# ---- 2.1 连续结局：分层回归 ----
for (oc in c("phq9_total", "gad7_total")) {
  dd <- p %>% select(all_of(c(oc, GOLD_VARS, ALPHA_VARS))) %>% drop_na()
  if (nrow(dd) < 40) next
  f0 <- as.formula(paste(oc, "~", paste(GOLD_VARS, collapse = " + ")))
  f1 <- as.formula(paste(oc, "~", paste(c(GOLD_VARS, ALPHA_VARS), collapse = " + ")))
  m0 <- lm(f0, data = dd); m1 <- lm(f1, data = dd)
  a <- anova(m0, m1)
  r2_0 <- summary(m0)$r.squared; r2_1 <- summary(m1)$r.squared
  cat(sprintf("  %s：金标准单独 R²=%.3f → 加入四 α 后 R²=%.3f（ΔR²=%+.3f，+%.0f%%）｜F(%d,%d)=%.2f，p=%s\n",
              oc, r2_0, r2_1, r2_1 - r2_0, (r2_1 - r2_0) / r2_0 * 100,
              a$Df[2], a$Res.Df[2], a$F[2], format.pval(a$`Pr(>F)`[2], digits = 3)))
  add_inc(oc, "连续（分层回归）", r2_0, r2_1,
          data.frame(n = nrow(dd), test_stat = round(a$F[2], 3),
                     p = format.pval(a$`Pr(>F)`[2], digits = 3),
                     AUC_gold_only = NA_real_, AUC_gold_plus_alpha = NA_real_,
                     delta_AUC = NA_real_, delong_p = NA_character_,
                     stringsAsFactors = FALSE))
}

# ---- 2.2 二分类结局：ΔAUC + DeLong ----
for (oc in c("phq9_pos", "gad7_pos")) {
  dd <- p %>% select(all_of(c(oc, GOLD_VARS, ALPHA_VARS))) %>% drop_na()
  if (nrow(dd) < 40 || length(unique(dd[[oc]])) < 2) next
  f0 <- as.formula(paste(oc, "~", paste(GOLD_VARS, collapse = " + ")))
  f1 <- as.formula(paste(oc, "~", paste(c(GOLD_VARS, ALPHA_VARS), collapse = " + ")))
  m0 <- suppressWarnings(glm(f0, data = dd, family = binomial))
  m1 <- suppressWarnings(glm(f1, data = dd, family = binomial))
  pr0 <- suppressWarnings(predict(m0, type = "response"))
  pr1 <- suppressWarnings(predict(m1, type = "response"))
  r0 <- pROC::roc(dd[[oc]], pr0, quiet = TRUE)
  r1 <- pROC::roc(dd[[oc]], pr1, quiet = TRUE)
  a0 <- as.numeric(pROC::auc(r0)); a1 <- as.numeric(pROC::auc(r1))
  dt <- tryCatch(pROC::roc.test(r0, r1, method = "delong"), error = function(e) NULL)
  pv <- if (is.null(dt)) NA_character_ else format.pval(dt$p.value, digits = 3)
  r2_0 <- 1 - m0$deviance / m0$null.deviance
  r2_1 <- 1 - m1$deviance / m1$null.deviance
  cat(sprintf("  %s：金标准单独 AUC=%.3f → 加入四 α 后 AUC=%.3f（ΔAUC=%+.3f）｜DeLong p=%s\n",
              oc, a0, a1, a1 - a0, pv))
  cat(sprintf("        伪 R²（McFadden 型）：%.3f → %.3f（ΔR²=%+.3f）\n", r2_0, r2_1, r2_1 - r2_0))
  add_inc(oc, "二分类（ΔAUC + DeLong）", r2_0, r2_1,
          data.frame(n = nrow(dd), test_stat = NA_real_, p = pv,
                     AUC_gold_only = round(a0, 4), AUC_gold_plus_alpha = round(a1, 4),
                     delta_AUC = round(a1 - a0, 4), delong_p = pv,
                     stringsAsFactors = FALSE))
}

if (length(irows)) save_tab(do.call(rbind, irows), "T22_incremental_validity")

cat("\n  ★ 判读要点：\n")
cat("    · ΔR²/ΔAUC 只描述同样本增量；不显著时不得宣称增量效度，显著时也不证明独立构念或因果机制\n")
cat("    · 表观 AUC 有乐观偏倚，CV 仍属内部验证；题量差异本身不构成等效或可替代证据\n")
cat("    · 收敛-区分判据须按上方实际比较结果报告；目标配对未系统更强时不得宣称该结构成立\n")
}
