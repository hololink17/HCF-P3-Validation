# =============================================================
# 10_analysis_00_common.R —— 分析层公共模块（两队列共用入口）
# -------------------------------------------------------------
# 作用：把 Prolific 与见数两队列读成同一套公开列名与下游统计分组。
# 只读取九个后台直出读数、外部变量与数值锚点，不重构内部计分。
# Prolific 每人一行；见数使用公开分析基底。
# =============================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(pROC); library(boot); library(ggplot2)
})

P3_SCRIPT_DIR <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L]), winslash = "/", mustWork = FALSE))
  } else {
    getwd()
  }
})
SELF_DIR <- P3_SCRIPT_DIR
PROJECT_ROOT <- normalizePath(Sys.getenv(
  "P3_PROJECT_ROOT", unset = file.path(SELF_DIR, "..", "..")
), winslash = "/", mustWork = FALSE)
BASE_DIR  <- PROJECT_ROOT
REPO_DIR  <- PROJECT_ROOT
R_SRC_DIR <- SELF_DIR
DERIVED_DIR <- normalizePath(Sys.getenv(
  "P3_DERIVED_DIR", unset = file.path(PROJECT_ROOT, "data")
), winslash = "/", mustWork = FALSE)
# 统一分析基底：两平台同构，由 build_analysis_base.py 生成。
# 公开管线只读取单一固定后台直出读数，不提供替代映射或替代参数开关。
base_path <- function(stem) {
  file.path(REPO_DIR, "code", "data", "derived", paste0(stem, ".csv"))
}
PROLIFIC_F <- base_path("p3_base_prolific")
JIANSHU_F  <- base_path("p3_base_jianshu")
GOLD_F     <- base_path("p3_base_gold")   # Prolific 完整版金标准子样本 n=142
DEDUP      <- TRUE          # TRUE = 每人一行（462）；FALSE = 沿用旧口径（471）

# ---------- 平台化运行（2026-09-29）----------
# P3_COHORT 决定「本次运行处理哪个队列 + 结果写到哪个目录」：
#   "prolific" → 只跑 Prolific 队列，结果写 BASE_DIR/数据分析/分析结果
#   "jianshu"  → 只跑见数队列，  结果写 BASE_DIR/见数/分析结果
#   "both"     → 两个队列同跑（缺省），结果写 PROJECT_ROOT/output
# 可用 P3_OUT_DIR 显式指定隔离输出目录；必须与 run_all 的清单/日志目录一致
# 设定方式（优先级从高到低）：
#   ① 调用方（run_all.R）已定义的全局变量 P3_COHORT
#   ② 环境变量 P3_COHORT
#   ③ 缺省 "both"
if (!exists("P3_COHORT", inherits = TRUE) || !nzchar(get0("P3_COHORT", inherits = TRUE))) {
  P3_COHORT <- Sys.getenv("P3_COHORT", "")
}
if (!nzchar(P3_COHORT)) P3_COHORT <- "both"
P3_COHORT <- tolower(trimws(P3_COHORT))

OUT_DIR <- Sys.getenv("P3_OUT_DIR", unset = switch(P3_COHORT,
  prolific = file.path(BASE_DIR, "数据分析", "分析结果"),
  jianshu  = file.path(BASE_DIR, "见数", "分析结果"),
  file.path(REPO_DIR, "output")))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
source(file.path(SELF_DIR, "00_shared_rules.R"))   # 唯一口径真值（与本文件同目录）

# ---------- 1. 统一基底读取（两平台同构，2026-10-26）----------
# ⚠ 架构变更：Prolific 与见数两平台现在读**同一套基底格式**
#   （代码/data/derived/p3_base_{prolific,jianshu}.csv，由 build_analysis_base.py 生成），
#   列名完全一致 → 下面的 read_base() 是**唯一**读取实现，两平台共用。
#   好处：口径不再可能分叉；官方量表计分（PHQ-9/GAD-7/MAIA-2/AAQ-II/DCS/BADS）
#   只在 Python 管线里实现一次，两平台必然同规则。
#
#   基底只含九个后台直出读数、外部公开量表结果、外部健康变量和数值锚点。
#   读取后仅额外派生书籍已公开的 H_trend = alpha1^2；不重构后台 H。
read_base <- function(f, cohort_label) {
  d <- read_csv(f, col_types = cols(.default = col_character()),
                show_col_types = FALSE) %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
  numc <- c("alpha1","alpha2","alpha3","alpha4",
            "carrier_a","carrier_b","carrier_c","carrier_d",
            "H_score",
            "phq9_total","gad7_total","gad7_item_count",
            "ibs_positive","ibs_symptom_count",
            "chronic_pain","heart_disease","hypertension","diabetes","hyperlipidemia",
            "bmi","obese","sleep_hour","smoke","alcohol_freq","alcohol_binary",
            "med_hypertension","med_diabetes","med_lipid","med_psych",
            "sleep_disorder_positive","sleep_symptom_count",
            "anchor1","anchor2","anchor3","anchor4",
            "maia_noticing","maia_not_distracting","maia_not_worrying","maia_three_dim",
            "aaq_total","dcs_total","bads_total",
            "duration_sec","depression","phq9_normal",
            "attention_pass","inattentive")
  d <- d %>% mutate(across(any_of(numc), ~ suppressWarnings(as.numeric(.x))))
  # 先在普通数据框环境检查可选列，避免 transmute 对不存在列产生初始化警告。
  d$flag_attention <- if ("inattentive" %in% names(d) && any(!is.na(d$inattentive)))
    d$inattentive else rep(NA_real_, nrow(d))
  d %>% transmute(
    cohort = cohort_label,
    pid    = if ("pid" %in% names(d)) .data$pid else .data$participant_token,
    data_source = .data$data_source,
    duration_sec = if ("duration_sec" %in% names(d)) .data$duration_sec else NA_real_,
    # 注意力标记：Prolific 用 invalid_*（0=未标记）；见数用 inattentive（1=误答）
    flag_attention = .data$flag_attention,
    flag_duration  = NA_real_,
    age_group = .data$age_group,
    gender    = .data$gender,
    education = .data$education,
    alpha1, alpha2, alpha3, alpha4,
    carrier_a, carrier_b, carrier_c, carrier_d,
    H_score,
    phq9_total, gad7_total,
    depression = .data$depression,          # 抑郁筛查阳性（PHQ-9 ≥ 10）
    ibs_positive, bmi,
    sleep_hour = .data$sleep_hour,
    sleep_symptom_count,
    sleep_disorder_positive,
    chronic_pain, heart_disease, hypertension, diabetes, hyperlipidemia,
    obese,
    smoke, alcohol_freq,
    med_hypertension, med_diabetes, med_lipid, med_psych,
    anchor1 = if ("anchor1" %in% names(d)) .data$anchor1 else NA_real_,
    anchor2 = if ("anchor2" %in% names(d)) .data$anchor2 else NA_real_,
    anchor3 = if ("anchor3" %in% names(d)) .data$anchor3 else NA_real_,
    anchor4 = if ("anchor4" %in% names(d)) .data$anchor4 else NA_real_,
    maia_three_dim, aaq_total, dcs_total, bads_total
  )
}

# ---------- 1a. Prolific（每人一行）----------
read_prolific <- function(dedup = DEDUP) {
  p <- read_base(PROLIFIC_F, "Prolific")
  if (dedup) {
    # 真实重复被试：同一 token 既参加完整版又参加精简版 → 优先保留 full
    # ⚠ 基底已按 submission_id 每人一行（462），此处保留接口以兼容 DEDUP=FALSE 口径
    p <- p %>% arrange(.data$pid, dplyr::desc(.data$data_source == "full")) %>%
      distinct(.data$pid, .keep_all = TRUE)
  }
  p
}

# ---------- 1b. 见数 ----------
read_jianshu <- function() {
  read_base(JIANSHU_F, "见数")
}

# ---------- 2b. 外部公开量表题目 ----------
# 公开 L2 不提供内部题目作答；缺少受控输入的分析应显式跳过。
read_prolific_items <- function(dedup = DEDUP) {
  p <- read_csv(PROLIFIC_F, col_types = cols(.default = col_character()),
                show_col_types = FALSE) %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
  numc <- c(paste0("alpha", 1:4), "H_score",
            paste0("phq_", 36:44), paste0("gad_", 45:51))
  p <- p %>% mutate(across(any_of(numc), ~ suppressWarnings(as.numeric(.x))))
  if (!"pid" %in% names(p) && "participant_token" %in% names(p)) p$pid <- p$participant_token
  if (!"data_source" %in% names(p)) p$data_source <- NA_character_
  if (dedup) {
    p <- p %>% arrange(.data$pid, dplyr::desc(.data$data_source == "full")) %>%
      distinct(.data$pid, .keep_all = TRUE)
  }
  p %>% transmute(
    cohort = "Prolific", pid = .data$pid, data_source = .data$data_source,
    alpha1, alpha2, alpha3, alpha4, H_score,
    across(any_of(c(paste0("phq_", 36:44), paste0("gad_", 45:51))))
  )
}

read_jianshu_items <- function() NULL

# ---------- 3. 载入 + 下游统计分组 ----------
# which: NULL（缺省）→ 按本次运行模式 P3_COHORT 决定读哪些队列
#        "both"      → 两队列都读
#        "prolific"  → 只读 Prolific（见数侧运行时不加载 Prolific 数据）
#        "jianshu"   → 只读见数
load_cohorts <- function(which = NULL) {
  if (is.null(which)) which <- P3_COHORT
  L <- list()
  if (which %in% c("both", "prolific"))
    L$prolific <- add_h_terms(add_alpha_typing(read_prolific()))
  if (which %in% c("both", "jianshu"))
    L$jianshu  <- add_h_terms(add_alpha_typing(read_jianshu()))
  L
}

# 跨队列对照专用：强制读入两个队列（无论本次运行模式）
load_both <- function() {
  list(prolific = add_h_terms(add_alpha_typing(read_prolific())),
       jianshu  = add_h_terms(add_alpha_typing(read_jianshu())))
}

# ---------- 4. 统计工具 ----------
# 单变量 AUC：按 glm 方向自动定向（等价于 R 的 roc(y, fitted(glm))）
auc_oriented <- function(y, x) {
  ok <- is.finite(y) & is.finite(x)
  y <- y[ok]; x <- x[ok]
  if (length(unique(y)) < 2) return(c(auc = NA, n = length(y), pos = sum(y)))
  a <- as.numeric(pROC::auc(pROC::roc(y, x, quiet = TRUE)))
  c(auc = max(a, 1 - a), n = length(y), pos = sum(y))
}

# AUC + 95% CI（DeLong）
auc_ci <- function(y, x, orient = TRUE) {
  ok <- is.finite(y) & is.finite(x); y <- y[ok]; x <- x[ok]
  r <- pROC::roc(y, if (orient) -x else x, quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(r))
  c(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3], n = length(y), pos = sum(y))
}

# 多变量模型 AUC（k 折交叉验证，手写实现，不依赖 caret —— caret 在本机 R 4.5.2 会段错误）
cv_auc <- function(df, preds, ycol = "ibs_positive", k = 5, seed = 42) {
  d <- df %>% select(all_of(c(ycol, preds))) %>% drop_na()
  if (nrow(d) < 30 || length(unique(d[[ycol]])) < 2) return(c(mean = NA, sd = NA, n = nrow(d)))
  f <- as.formula(paste(ycol, "~", paste(preds, collapse = "+")))
  set.seed(seed)
  folds <- sample(rep(1:k, length.out = nrow(d)))
  aucs <- c()
  for (i in 1:k) {
    tr <- d[folds != i, ]; te <- d[folds == i, ]
    if (length(unique(tr[[ycol]])) < 2 || length(unique(te[[ycol]])) < 2) next
    m <- suppressWarnings(glm(f, data = tr, family = binomial))
    pr <- suppressWarnings(predict(m, newdata = te, type = "response"))
    aucs <- c(aucs, as.numeric(pROC::auc(pROC::roc(te[[ycol]], pr, quiet = TRUE))))
  }
  c(mean = mean(aucs), sd = sd(aucs), n = nrow(d), k = length(aucs))
}

# 两队列 AUC 差异的 bootstrap CI（DeLong 的稳健替代实现）
auc_diff_boot <- function(y1, x1, y2, x2, B = 3000, seed = 7) {
  ok1 <- is.finite(y1) & is.finite(x1); y1 <- y1[ok1]; x1 <- x1[ok1]
  ok2 <- is.finite(y2) & is.finite(x2); y2 <- y2[ok2]; x2 <- x2[ok2]
  observed <- as.numeric(pROC::auc(pROC::roc(y1, -x1, quiet = TRUE, direction = "<"))) -
    as.numeric(pROC::auc(pROC::roc(y2, -x2, quiet = TRUE, direction = "<")))
  set.seed(seed)
  one <- function() {
    i1 <- sample(seq_along(y1), replace = TRUE); i2 <- sample(seq_along(y2), replace = TRUE)
    if (length(unique(y1[i1])) < 2 || length(unique(y2[i2])) < 2) return(NA_real_)
    a1 <- tryCatch(as.numeric(pROC::auc(pROC::roc(y1[i1], -x1[i1], quiet = TRUE, direction = "<"))), error = function(e) NA_real_)
    a2 <- tryCatch(as.numeric(pROC::auc(pROC::roc(y2[i2], -x2[i2], quiet = TRUE, direction = "<"))), error = function(e) NA_real_)
    a1 - a2
  }
  d <- vapply(seq_len(B), function(i) one(), numeric(1))
  d <- d[is.finite(d)]
  if (length(d) < 100) return(c(delta = NA_real_, lo = NA_real_, hi = NA_real_, nboot = length(d)))
  c(delta = observed, lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975)), nboot = length(d))
}

# 效应量：Cohen's d
cohens_d <- function(x, g) {
  a <- x[g == 1]; b <- x[g == 0]
  if (length(a) < 3 || length(b) < 3) return(NA_real_)
  s <- sqrt(((length(a)-1)*var(a) + (length(b)-1)*var(b)) / (length(a)+length(b)-2))
  (mean(a) - mean(b)) / s
}

# 三行式结果落盘
save_tab <- function(x, name) {
  f <- file.path(OUT_DIR, paste0(name, ".csv"))
  readr::write_excel_csv(as.data.frame(x), f)
  cat("  → 已写出:", basename(f), "\n")
  invisible(f)
}

cat("=== 10_analysis_00_common.R 载入完成 ===\n")
cat("  运行模式 P3_COHORT:", P3_COHORT,
    switch(P3_COHORT, prolific = "（仅 Prolific 队列）", jianshu = "（仅见数队列）",
           "（两队列同跑）"), "\n")
cat("  结果输出目录:", OUT_DIR, "\n")
cat("  Prolific 口径:", ifelse(DEDUP, "每人一行", "含重复行"), "\n")
cat("  下游统计分组规则源: 00_shared_rules.R\n")
cat("  版本:", VERSION_TAG, "（两队列一致）\n")
