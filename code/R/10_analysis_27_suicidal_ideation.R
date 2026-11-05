# =============================================================
# 10_analysis_27_suicidal_ideation.R —— H8：α₁ 与 PHQ-9 第 9 条（自杀意念条目）
# -------------------------------------------------------------
# 背景（2026-10-04）：
#   此前判定 H8「结构性不可做、量表未采集自杀意念条目」—— **该判定已作废**。
#   错因：在问卷 docx 内检索关键词「自杀」零命中，据此下了**否定性结论**；
#   而问卷原文第 9 条写作「想到自己最好去死或者自残」，不含「自杀」二字。
#   🔴 铁律：否定性结论（「没有 X」）不得由关键词检索得出，必须穷举或读原文结构。
#
# 三重编号印证（两队列都采集了该条目）：
#   ① 数据字典 `docs/分析用数据字典_v1_20260927.md` §5.3：
#        PHQ-9 → Prolific `phq_36`…`phq_44` ｜ 见数 `phq_25`…`phq_33`
#   ② 见数 `代码2_构建分析基底.R` L106/L109：
#        phq_items <- paste0("Q", 25:33); phq9_N <- likert(Q(N+24))
#        → **phq9_9 ← Q33**
#   ③ 见数问卷原文 `调查问卷中文60题2.0.docx`：Q33 =「想到自己最好去死或者自残」
#      Prolific 侧 PHQ-9 第 9 条原文（`题目/PHQ-9.docx` 表）：
#        "Thoughts that you would be better off dead or of hurting yourself in some way"
#
# 本模块四节：
#   27.0 条目位置核验（9 题剖面）—— 防「取错列」的守卫，不作为结果报告
#   27.1 H8 主检验：Spearman r(α₁, item9)，双队列
#   27.2 二分类对照：item9 > 0 vs = 0 的 α₁ 差异（Cohen's d ＋ logistic OR）
#   27.3 α₁ 特异性：四个 α 各自与 item9 的关联，看 α₁ 是否最强
#   27.4 独立于抑郁整体：控制 PHQ-9 其余 8 题后的偏相关
#   27.5 对照：α₁ vs PHQ-9 总分
#
# 🔴 写作口径（三条必须同现，否则属过度解读）：
#   ① 主检验支持（双队列共同显著）——本项目罕见
#   ② **非特异**：控制其余 8 题后偏相关不显著 → 该关联可由整体抑郁解释
#   ③ α₁ **不特异**：四个 α 效应量相当（Prolific 最强为 α₄）→ 与 H2 区分效度失败同源
#   故论文写法＝「H8 在零阶水平得到支持，但不具备条目特异性」，
#   **禁写**「工具可识别自杀风险」「α₁ 可预测自杀意念」等临床可用性主张。
# =============================================================
# --- Windows 中文路径守卫 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

ALPHA_LAB <- c(a1 = "α₁", a2 = "α₂", a3 = "α₃", a4 = "α₄")
ITEM_LAB  <- paste0("PHQ-9 第 ", 1:9, " 条")

cat("\n=== H8：α₁ 与 PHQ-9 第 9 条（自杀意念条目）===\n")
cat(sprintf("队列 = %s\n", P3_COHORT))

# -------------------------------------------------------------
# 取数层
# 🔴 关键：去重口径必须与 10_analysis_00_common.R 的 read_prolific() **逐字一致**
#    （arrange(participant_token, desc(data_source == "full"))）。
#    prolific_clean.csv 有 471 行 / 462 个唯一 token；9 行重复中
#    **4 行是 full+short 混合**（须保留 full）、5 行是 short+short。
#    若改用「按文件出现顺序保留首行」，会在这 4 人上取到 short 版数值 → 与主口径不符。
# -------------------------------------------------------------
read_items_prolific <- function() {
  p <- base_path("p3_base_prolific")
  d <- readr::read_csv(p, col_types = readr::cols(.default = readr::col_character()),
                       show_col_types = FALSE) %>%
    dplyr::mutate(dplyr::across(dplyr::where(is.character), ~ dplyr::na_if(.x, "")))
  numc <- c(paste0("alpha", 1:4), paste0("phq_", 36:44))
  d <- d %>% dplyr::mutate(dplyr::across(dplyr::any_of(numc),
                                         ~ suppressWarnings(as.numeric(.x))))
  req <- c("participant_token", "data_source", paste0("alpha", 1:4), paste0("phq_", 36:44))
  if (!all(req %in% names(d))) return(NULL)
  if (DEDUP) {
    d <- d %>% dplyr::arrange(.data$participant_token,
                              dplyr::desc(.data$data_source == "full")) %>%
      dplyr::distinct(.data$participant_token, .keep_all = TRUE)
  }
  d %>% dplyr::transmute(
    cohort = "Prolific", pid = .data$participant_token,
    a1 = .data$alpha1, a2 = .data$alpha2, a3 = .data$alpha3, a4 = .data$alpha4,
    it1 = .data$phq_36, it2 = .data$phq_37, it3 = .data$phq_38,
    it4 = .data$phq_39, it5 = .data$phq_40, it6 = .data$phq_41,
    it7 = .data$phq_42, it8 = .data$phq_43, it9 = .data$phq_44
  )
}

read_items_jianshu <- function() {
  j <- readr::read_csv(JIANSHU_F, col_types = readr::cols(.default = readr::col_character()),
                       show_col_types = FALSE)
  # ⚠ 2026-10-26 统一基底口径：α 列名为 alpha1..alpha4、PHQ 条目列为 phq_36..phq_44
  #   （旧 p3_style_analysis_base.csv 用 scores_alpha* / phq9_1..9）
  numc <- c(paste0("alpha", 1:4), paste0("phq_", 36:44))
  j <- j %>% dplyr::mutate(dplyr::across(dplyr::any_of(numc),
                                         ~ suppressWarnings(as.numeric(.x))))
  id_col <- if ("pid" %in% names(j)) "pid" else "participant_token"
  req <- c(id_col, paste0("alpha", 1:4), paste0("phq_", 36:44))
  if (!all(req %in% names(j))) return(NULL)
  j %>% dplyr::transmute(
    cohort = "见数", pid = .data[[id_col]],
    a1 = .data$alpha1, a2 = .data$alpha2,
    a3 = .data$alpha3, a4 = .data$alpha4,
    it1 = .data$phq_36, it2 = .data$phq_37, it3 = .data$phq_38,
    it4 = .data$phq_39, it5 = .data$phq_40, it6 = .data$phq_41,
    it7 = .data$phq_42, it8 = .data$phq_43, it9 = .data$phq_44
  )
}

# -------------------------------------------------------------
# 工具函数
# -------------------------------------------------------------
# Spearman：cor.test 的精确 p 在存在并列秩时不可用 → exact = FALSE（AS89 近似）
sp_test <- function(x, y) {
  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  c(r = unname(ct$estimate), p = ct$p.value)
}

# Cohen's d（合并 SD）
cohen_d <- function(g1, g0) {
  n1 <- length(g1); n0 <- length(g0)
  if (n1 < 2 || n0 < 2) return(NA_real_)
  sp <- sqrt(((n1 - 1) * stats::var(g1) + (n0 - 1) * stats::var(g0)) / (n1 + n0 - 2))
  if (!is.finite(sp) || sp == 0) return(NA_real_)
  (mean(g1) - mean(g0)) / sp
}

# 偏相关：控制 Z（data.frame，每列一个协变量）后 x 与 y 的相关
# df = n − k − 2（k = 协变量个数）—— 这是标准偏相关检验自由度
pcor_manual <- function(x, y, Z) {
  n <- length(x); k <- ncol(Z)
  rx <- stats::residuals(stats::lm(x ~ ., data = Z))
  ry <- stats::residuals(stats::lm(y ~ ., data = Z))
  r <- stats::cor(rx, ry)
  df <- n - k - 2
  if (df < 1 || !is.finite(r) || abs(r) >= 1) return(c(r = r, t = NA_real_, df = df, p = NA_real_))
  t <- r * sqrt(df / (1 - r^2))
  c(r = r, t = t, df = df, p = 2 * stats::pt(-abs(t), df))
}

# -------------------------------------------------------------
# 主体
# -------------------------------------------------------------
run_cohort <- function(d, nm) {
  d <- d %>% dplyr::filter(dplyr::if_all(dplyr::all_of(c(paste0("it", 1:9),
                                                          paste0("a", 1:4))),
                                         ~ !is.na(.x)))
  n <- nrow(d)
  if (n < 40) {
    cat(sprintf("  [%s] n=%d 过小，跳过\n", nm, n))
    return(NULL)
  }
  cat(sprintf("\n  [%s] 有效 n = %d\n", nm, n))

  it9 <- d$it9

  # ---------- 27.0 条目位置核验（守卫）----------
  prof <- data.frame(
    cohort = nm, item = ITEM_LAB,
    mean = round(sapply(d[paste0("it", 1:9)], mean), 3),
    sd = round(sapply(d[paste0("it", 1:9)], stats::sd), 3),
    pct_gt0 = round(sapply(d[paste0("it", 1:9)], function(v) mean(v > 0) * 100), 1),
    row.names = NULL
  )
  cat("  --- 27.0 PHQ-9 九题剖面（守卫：第 9 条应为低分高偏态项）---\n")
  print(prof)
  cat(sprintf("  第 9 条阳性率 = %.1f%%（>0 者 %d 人）\n",
              mean(it9 > 0) * 100, sum(it9 > 0)))

  # ---------- 27.1 H8 主检验 ----------
  s9 <- sp_test(d$a1, it9)
  cat(sprintf("  27.1 H8 主检验：Spearman r(α₁, item9) = %+.4f  p = %.4g  n = %d\n",
              s9["r"], s9["p"], n))

  # ---------- 27.2 二分类对照 ----------
  g1 <- d$a1[it9 > 0]; g0 <- d$a1[it9 == 0]
  dd <- cohen_d(g1, g0)
  tt <- tryCatch(suppressWarnings(stats::t.test(g1, g0)), error = function(e) NULL)
  ww <- tryCatch(suppressWarnings(stats::wilcox.test(g1, g0)), error = function(e) NULL)
  ord <- tryCatch(suppressWarnings(stats::glm(it9 > 0 ~ d$a1, family = stats::binomial)),
                  error = function(e) NULL)
  or_crude <- if (is.null(ord)) NA_real_ else exp(unname(stats::coef(ord)[2]))
  or_ci <- if (is.null(ord)) c(NA_real_, NA_real_) else
    suppressMessages(exp(stats::confint(ord)[2, ]))
  or_sd <- if (is.null(ord)) NA_real_ else exp(unname(stats::coef(ord)[2]) * stats::sd(d$a1))
  cat(sprintf("  27.2 二分类：α₁ 均值 >0 组 = %.3f ｜ =0 组 = %.3f ｜ 差 = %+.3f ｜ Cohen's d = %+.3f\n",
              mean(g1), mean(g0), mean(g1) - mean(g0), dd))
  if (!is.null(tt)) cat(sprintf("       t 检验 p = %.4g ｜ Wilcoxon p = %.4g\n", tt$p.value, ww$p.value))
  cat(sprintf("       粗 OR（每 +1 α₁）= %.4f（95%%CI %.4f–%.4f）｜ 每 +1 SD = %.4f\n",
              or_crude, or_ci[1], or_ci[2], or_sd))

  # ---------- 27.3 α₁ 特异性：四个 α 对照 ----------
  cmp <- data.frame(
    cohort = nm,
    alpha = ALPHA_LAB,
    r = NA_real_, p = NA_real_,
    row.names = NULL
  )
  for (a in names(ALPHA_LAB)) {
    s <- sp_test(d[[a]], it9)
    cmp$r[cmp$alpha == ALPHA_LAB[a]] <- unname(s["r"])
    cmp$p[cmp$alpha == ALPHA_LAB[a]] <- unname(s["p"])
  }
  cmp$r <- round(cmp$r, 4)
  cmp$p <- signif(cmp$p, 3)
  cmp$is_h8 <- ifelse(cmp$alpha == "α₁", "← H8", "")
  cat("  27.3 四个 α 分别与 item9 的 Spearman（检验 α₁ 是否特异）---\n")
  print(cmp)
  strongest <- cmp$alpha[which.min(cmp$r)]
  cat(sprintf("       最强负向 = %s（r = %+.4f）→ %s\n", strongest,
              min(cmp$r),
              ifelse(strongest == "α₁", "α₁ 最强", "🔴 α₁ **不是**最强 → 无条目特异性")))

  # ---------- 27.4 独立于抑郁整体：偏相关 ----------
  Z <- as.data.frame(d[paste0("it", 1:8)])          # 控制其余 8 题
  pc <- pcor_manual(d$a1, it9, Z)
  cat(sprintf("  27.4 偏相关 r(α₁, item9 | 其余 8 题) = %+.4f  t = %+.3f  df = %d  p = %.4g\n",
              pc["r"], pc["t"], pc["df"], pc["p"]))
  cat(sprintf("       → %s\n",
              ifelse(pc["p"] < 0.05,
                     "✅ 独立于抑郁整体",
                     "⚠️ 不显著 → 该关联可由整体抑郁解释，**无条目特异性贡献**")))

  # ---------- 27.5 对照：α₁ vs PHQ-9 总分 ----------
  tot <- rowSums(d[paste0("it", 1:9)])
  st <- sp_test(d$a1, tot)
  cat(sprintf("  27.5 对照：Spearman r(α₁, PHQ-9 总分) = %+.4f  p = %.4g\n", st["r"], st["p"]))

  list(
    main = data.frame(
      cohort = nm, n = n,
      n_item9_pos = sum(it9 > 0),
      pct_item9_pos = round(mean(it9 > 0) * 100, 1),
      spearman_r = round(s9["r"], 4),
      p = signif(s9["p"], 3),
      a1_mean_pos = round(mean(g1), 3),
      a1_mean_neg = round(mean(g0), 3),
      a1_diff = round(mean(g1) - mean(g0), 3),
      cohen_d = round(dd, 3),
      ttest_p = if (!is.null(tt)) signif(tt$p.value, 3) else NA_real_,
      wilcox_p = if (!is.null(ww)) signif(ww$p.value, 3) else NA_real_,
      OR_crude = round(or_crude, 4),
      OR_lo = round(or_ci[1], 4), OR_hi = round(or_ci[2], 4),
      OR_per_SD = round(or_sd, 4),
      partial_r = round(pc["r"], 4),
      partial_df = pc["df"],
      partial_p = signif(pc["p"], 3),
      r_phq_total = round(st["r"], 4),
      p_phq_total = signif(st["p"], 3)
    ),
    cmp = cmp,
    prof = prof
  )
}

res <- list()
missing_l3 <- FALSE
if (P3_COHORT %in% c("prolific", "both")) {
  d <- read_items_prolific()
  if (is.null(d)) {
    missing_l3 <- TRUE
    cat("  Prolific：公开 L2 不含 PHQ-9 原始题项，跳过；本分析需要受控 L3。\n")
  } else {
    r <- tryCatch(run_cohort(d, "prolific"), error = function(e) {
      cat("  ⚠ Prolific 侧失败：", conditionMessage(e), "\n")
      .P3_MODULE_STATUS <<- "failed"
      .P3_MODULE_REASON <<- conditionMessage(e)
      NULL })
    if (!is.null(r)) res[["prolific"]] <- r
  }
}
if (P3_COHORT %in% c("jianshu", "both")) {
  d <- read_items_jianshu()
  if (is.null(d)) {
    missing_l3 <- TRUE
    cat("  见数：公开 L2 不含 PHQ-9 原始题项，跳过；本分析需要受控 L3。\n")
  } else {
    r <- tryCatch(run_cohort(d, "jianshu"), error = function(e) {
      cat("  ⚠ 见数侧失败：", conditionMessage(e), "\n")
      .P3_MODULE_STATUS <<- "failed"
      .P3_MODULE_REASON <<- conditionMessage(e)
      NULL })
    if (!is.null(r)) res[["jianshu"]] <- r
  }
}

if (!length(res)) {
  if (missing_l3 && (!exists(".P3_MODULE_STATUS") || !identical(.P3_MODULE_STATUS, "failed"))) {
    .P3_MODULE_STATUS <- "skipped_controlled_l3_required"
    .P3_MODULE_REASON <- "需要受控 L3 的 PHQ-9 原始题项"
  }
  cat("  ⚠ 无可用队列，退出\n")
} else {
  if (!identical(.P3_MODULE_STATUS, "failed") && length(res) < if (identical(P3_COHORT, "both")) 2L else 1L) {
    .P3_MODULE_STATUS <- "partial_controlled_l3"
    .P3_MODULE_REASON <- sprintf("仅完成%s；其余队列缺少受控题项或未成功", paste(names(res), collapse = "+"))
  }
  main_all <- do.call(rbind, lapply(res, `[[`, "main"))
  write.csv(main_all, file.path(OUT_DIR, "T27_h8_suicidal_ideation.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n[写出] T27_h8_suicidal_ideation.csv\n")
  print(main_all)

  cmp_all <- do.call(rbind, lapply(res, `[[`, "cmp"))
  write.csv(cmp_all, file.path(OUT_DIR, "T27_h8_alpha_compare.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("[写出] T27_h8_alpha_compare.csv\n")

  prof_all <- do.call(rbind, lapply(res, `[[`, "prof"))
  write.csv(prof_all, file.path(OUT_DIR, "T27_phq9_item_profile.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("[写出] T27_phq9_item_profile.csv\n")

  # ---------- 结论口径 ----------
  cat("\n--- 结论口径 ---\n")
  for (nm in names(res)) {
    m <- res[[nm]]$main
    cat(sprintf(
      "%s：n=%d｜item9 阳性 %.1f%%｜r(α₁,item9)=%+.4f (p=%.3g)｜偏相关 %+.4f (p=%.3g)｜r(α₁,PHQ总分)=%+.4f\n",
      nm, m$n, m$pct_item9_pos, m$spearman_r, m$p, m$partial_r, m$partial_p, m$r_phq_total))
  }
  cat("\n  📌 写作三条同现：① 零阶支持（双队列）② 控制后非特异 ③ 四个 α 效应量相当（α₁ 不特异）\n")
  cat("  🔴 禁写「工具可识别自杀风险」「α₁ 可预测自杀意念」等临床可用性主张。\n")
  cat("  🔴 该条目属敏感议题，投稿时按目标期刊要求考虑在 Limitation 说明数据性质与安全预案。\n")
}
