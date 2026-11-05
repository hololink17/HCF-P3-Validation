# =============================================================
# 10_analysis_17_screening_equivalence.R —— 20 题 vs 85 题金标准的筛查模型对照
# 注意：本模块执行差异检验与内部 CV，不执行 TOST 等效性检验；历史文件/表名不代表等效结论。
# -------------------------------------------------------------
# 对应蓝图 R4：HCF（20 题）与金标准组合（85 题）在抑郁筛查上是否等效
# 仅 Prolific 完整版可做（n≈142，需四项金标准同时在测）
#
# 金标准组合构成（共 85 题，可核验）：
#   MAIA-2  37 题 → 3 个维度分（Noticing / Not-Distracting / Not-Worrying）
#   AAQ-II   7 题 → aaq_total
#   DCS     16 题 → dcs_total
#   BADS    25 题 → bads_total
#
# 模型（结局：PHQ-9 ≥ 10，切点先验指定）：
#   A  HCF：四 α（20 题）
#   B  金标准：上述 6 个已计分维度（85 题）
#   C  HCF + 金标准（105 题）
#
# 统计口径（面向审稿人）：
#   ① 主检验用 DeLong 检验比较相关 ROC 曲线（同一批样本、同一结局）
#   ② 同时给出 5 折交叉验证 AUC 作为内部验证，可减轻但不消除表观乐观偏倚
#   ③ A、B 均嵌套于 C；嵌套模型另附似然比检验（anova）作为补充
#   ④ DeLong 未显著不等于等效；本模块不主张等效、替代或超越金标准
#
# 输出：T14_screening_equivalence.csv
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
num <- c(ALPHA_AC, "phq9_total", "gad7_total",
         "maia_noticing", "maia_not_distracting", "maia_not_worrying",
         "maia_three_dim", "aaq_total", "dcs_total", "bads_total")

g <- read_csv(GOLD_F, col_types = cols(.default = col_character()), show_col_types = FALSE) %>%
  mutate(across(any_of(num), ~ suppressWarnings(as.numeric(.x)))) %>%
  filter(!is.na(maia_three_dim), !is.na(aaq_total), !is.na(dcs_total), !is.na(bads_total))

n_raw <- nrow(g)
if (DEDUP) {
  g <- g %>% arrange(participant_token, desc(data_source == "full")) %>%
    distinct(participant_token, .keep_all = TRUE)
}
cat(sprintf("\n=== 筛查模型对照（非等效性检验；完整版 %d 行 → 去重后 %d 行）===\n", n_raw, nrow(g)))

g$dep <- as.numeric(g$phq9_total >= 10)
cat(sprintf("  抑郁阳性 PHQ-9≥10：%d / %d（%.1f%%）\n",
            sum(g$dep, na.rm = TRUE), nrow(g), 100 * mean(g$dep, na.rm = TRUE)))

HCF_P  <- ALPHA_AC
GOLD_P <- c("maia_noticing", "maia_not_distracting", "maia_not_worrying",
            "aaq_total", "dcs_total", "bads_total")
stopifnot(all(c(HCF_P, GOLD_P, "dep") %in% names(g)))
g_m <- g[complete.cases(g[, c("dep", HCF_P, GOLD_P)]), ]
cat(sprintf("  建模完整案例 n=%d\n", nrow(g_m)))

fit_glm <- function(preds) {
  f <- as.formula(paste("dep ~", paste(preds, collapse = " + ")))
  suppressWarnings(glm(f, data = g_m, family = binomial))
}
m_A <- fit_glm(HCF_P)
m_B <- fit_glm(GOLD_P)
m_C <- fit_glm(c(HCF_P, GOLD_P))

roc_of <- function(m) pROC::roc(g_m$dep, suppressWarnings(predict(m, type = "response")), quiet = TRUE)
rA <- roc_of(m_A); rB <- roc_of(m_B); rC <- roc_of(m_C)

auc_row <- function(r, tag, nq) {
  ci <- as.numeric(pROC::ci.auc(r))
  data.frame(model = tag, n_items = nq,
             auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3],
             n = length(r$response), pos = sum(r$response), row.names = NULL)
}
res <- do.call(rbind, list(
  auc_row(rA, "A｜HCF 四 α", 20),
  auc_row(rB, "B｜金标准组合（MAIA-2 + AAQ-II + DCS + BADS）", 85),
  auc_row(rC, "C｜HCF + 金标准", 105)))

cat("\n  ── AUC ──\n"); print(res, row.names = FALSE)

# ---------- DeLong 检验 ----------
dl <- function(r1, r2, lab) {
  t <- pROC::roc.test(r1, r2, method = "delong")
  data.frame(comparison = lab,
             auc_1 = as.numeric(pROC::auc(r1)), auc_2 = as.numeric(pROC::auc(r2)),
             delta_auc = as.numeric(pROC::auc(r1)) - as.numeric(pROC::auc(r2)),
             Z = unname(t$statistic), p = t$p.value, test = "DeLong", row.names = NULL)
}
dl_tab <- rbind(
  dl(rA, rB, "A（HCF 20 题） vs B（金标准 85 题）"),
  dl(rC, rA, "C（HCF+金标准 105 题） vs A（HCF 20 题）"),
  dl(rC, rB, "C（HCF+金标准 105 题） vs B（金标准 85 题）"))

cat("\n  ── DeLong 检验 ──\n"); print(dl_tab, row.names = FALSE)

# ---------- 补充：嵌套模型的似然比检验（C vs A）----------
lr <- anova(m_A, m_C, test = "Chisq")
lr_row <- data.frame(
  comparison = "C vs A（嵌套，似然比检验）",
  delta_chi2 = lr$Deviance[2], delta_df = lr$Df[2],
  p = lr$`Pr(>Chi)`[2], row.names = NULL)
cat("\n  ── 嵌套模型似然比检验（DeLong 的补充）──\n"); print(lr_row, row.names = FALSE)

# ---------- 内部验证：5 折交叉验证 AUC ----------
cv_A <- cv_auc(transform(g_m, y = g_m$dep), HCF_P,  ycol = "y", k = 5)
cv_B <- cv_auc(transform(g_m, y = g_m$dep), GOLD_P, ycol = "y", k = 5)
cv_C <- cv_auc(transform(g_m, y = g_m$dep), c(HCF_P, GOLD_P), ycol = "y", k = 5)
cv_tab <- data.frame(
  model = c("A｜HCF 四 α（20 题）", "B｜金标准组合（85 题）", "C｜HCF + 金标准（105 题）"),
  cv_auc = c(cv_A["mean"], cv_B["mean"], cv_C["mean"]),
  cv_sd  = c(cv_A["sd"],   cv_B["sd"],   cv_C["sd"]), row.names = NULL)
cat("\n  ── 5 折交叉验证 AUC（内部验证；不是外部验证）──\n"); print(cv_tab, row.names = FALSE)

save_tab(res, "T14_screening_equivalence")
save_tab(dl_tab, "T14b_screening_delong")
save_tab(cv_tab, "T14c_screening_cv")
save_tab(lr_row, "T14d_screening_lrt")

cat("\n  ⚠️ 报告纪律：\n")
cat("     · 结局切点（PHQ-9≥10）先验指定；A/B/C 三模型预测变量集合先验指定\n")
cat("     · 表观 AUC 存在乐观偏倚；5 折 CV 是内部验证，可减轻但不消除偏倚，且不是外部验证\n")
cat("     · A、B 均嵌套于 C，故补嵌套似然比检验；DeLong 未显著不等于等效\n")
cat("     · 本模块不执行等效性检验，不主张等效、替代或超越金标准；题量差异须与准确度并列报告\n")
cat("     · 金标准仅完整版在测（问卷版本差异非随机分组），须在 Limitations 交代\n")

}
