# =============================================================
# 10_analysis_26_item_structure_efa.R —— R3 【零成本补强】20 题题目结构：EFA（无预设归属）
# -------------------------------------------------------------
# 为什么做 EFA 而不做 CFA（H3）：
#   ① CFA 必须预设「哪几题归哪个 α」——该映射属**专有资产**（内部映射表，不随仓库分发），
#      写进代码即触碰保密边界。EFA **不需要任何预设归属**，让数据自己说话。
#   ② H3 的 CFA 是**反身性（reflective）**测量模型的检验；而 HCF 四 α 在设计上是
#      **形成性／计算性指数（formative/computational index）**。以反身性模型检验形成性指数，
#      属**模型误配**，其失败不构成对工具的否定。
#   ③ 因此本模块的定位是：**如实报告题目层面的数据结构**，并据此**支持 formative 立场**，
#      而不是"试了确认性模型但失败"。
#
# 本模块输出（不含任何映射信息）：
#   26.1 平行分析 + 碎石图判据：保留几个因子（Horn 平行分析，比 Kaiser 特征值>1 稳健）
#   26.2 EFA 载荷矩阵（1–6 因子，oblimin 斜交旋转）→ 每题最大载荷归属
#   26.3 单因子性检验：首因子解释比、Cronbach α（20 题总量表）、平均题间相关
#   26.4 结论口径：四因子结构是否得到支持
#
# 🔴 保密纪律：本模块只在受控输入可用时读取题目级作答；公开 L2 不含这些字段，
#   因此公开运行不包含内部题目级分析，也不推断任何题目归属。
# =============================================================
# --- Windows 中文路径守卫 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

suppressPackageStartupMessages({
  library(psych)      # fa.parallel / principal / KMO
  library(ggplot2)
})

ITEM_COLS <- paste0("hcf_", 5:24)   # 20 题

cat("\n=== 20 题题目结构：EFA（无预设归属）===\n")
cat(sprintf("队列 = %s ｜ 共 %d 题\n", P3_COHORT, length(ITEM_COLS)))

# ---------- 取数 ----------
# 🔴 只做 Prolific 侧：见数侧的题目级数据为**选项文本**（cleaned_60q_merged 的 Q1–Q60），
#    公开 L2 不含内部题目作答；缺少受控输入时跳过该分支。
#    故本模块的题目级结构分析**仅限 Prolific**（分值口径，采集期两批标签下数值一致）。
#    见数侧的结构证据由模块 24（α 层面等价性）承担。
read_items_prolific <- function() {
  p <- base_path("p3_base_prolific")
  d <- readr::read_csv(p, col_types = readr::cols(.default = readr::col_character()),
                       show_col_types = FALSE)
  d %>% dplyr::mutate(across(where(is.character), ~ dplyr::na_if(.x, "")))
}

# 🔴 2026-10-04 修正：原写法用 quit(save="no") 结束进程 —— 本模块被 10_analysis_run_all.R
#    以 source() 在**同一进程**内调用，quit() 会**终止整个 run_all**
#    （见数模式下模块 27 及其后全部不执行，且 sink() 不关闭）。
#    改为条件包裹：跳过时脚本正常结束并交回 run_all。
SKIP_JIANSHU <- (P3_COHORT == "jianshu")
if (SKIP_JIANSHU) {
  cat("\n  ⚠ 本模块仅对 Prolific 队列执行；其他受控题目级输入不在此处读取。\n")
  cat("     见数侧结构证据见模块 24（α 层面等价性）。跳过。\n")
} else {

d <- read_items_prolific()
# 口径：每人一行（DEDUP=TRUE 时与主口径一致）——prolific_clean 已是清洗层，
# 🔴 去重必须与 10_analysis_00_common.R 的 read_prolific() **逐字一致**：
#    arrange(participant_token, desc(data_source == "full"))。
#    471 行含 9 行重复，其中 **4 行是 full+short 混合**，须保留 full 版；
#    若只按 token 排序去重（丢失 desc(full) 次序），这 4 人可能取到 short 版数值，
#    与主口径（462）不符。（2026-10-04 修正）
if (DEDUP && "participant_token" %in% names(d)) {
  d <- d %>% dplyr::arrange(.data$participant_token,
                            dplyr::desc(.data$data_source == "full")) %>%
    dplyr::distinct(.data$participant_token, .keep_all = TRUE)
}
cat(sprintf("去重后 n = %d\n", nrow(d)))

run_efa <- function(d, cohort_label, n_for_label) {
  X <- d %>% select(all_of(ITEM_COLS)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.x)))) %>%
    as.data.frame()
  X <- X[stats::complete.cases(X), , drop = FALSE]
  n <- nrow(X)
  if (n < 50) {
    cat(sprintf("  [%s] n=%d 过小，跳过\n", cohort_label, n))
    return(NULL)
  }
  cat(sprintf("  [%s] 有效 n = %d\n", cohort_label, n))

  # ---------- 26.1 平行分析：保留因子数 ----------
  # 🔴 兼容性：本机 psych 版本的 fa.parallel 不接受 quiet 参数，传了会整体报错
  #    （曾被 tryCatch 吞成 NULL → parallel_factors 恒为 NA）。故分两级调用。
  pa <- tryCatch(
    fa.parallel(X, fa = "fa", fm = "ml", n.iter = 100, plot = FALSE),
    error = function(e) tryCatch(
      fa.parallel(X, fa = "fa", fm = "ml", n.iter = 100, plot = FALSE, quiet = TRUE),
      error = function(e2) NULL
    )
  )
  # 返回字段随 psych 版本不同：$nfact（本机存在）/ $nfact_parallel / $nfactor → 逐级回退
  pick_nfact <- function(x) {
    if (is.null(x)) return(NA_integer_)
    for (k in c("nfact", "nfact_parallel", "nfactor")) {
      v <- x[[k]]
      if (!is.null(v) && length(v) >= 1 && !all(is.na(v))) return(as.integer(v[1]))
    }
    NA_integer_
  }
  nfact_pa <- pick_nfact(pa)
  # 兜底：用平行分析的「真实特征值 vs 随机特征值」曲线自行计数
  if (is.na(nfact_pa) && !is.null(pa)) {
    pcv <- tryCatch(as.numeric(pa$pc.values), error = function(e) NULL)
    fav <- tryCatch(as.numeric(pa$fa.values), error = function(e) NULL)
    if (!is.null(pcv) && !is.null(fav) && length(pcv) == length(fav)) {
      nfact_pa <- as.integer(sum(pcv > fav))
    }
  }
  # 特征值（相关矩阵）
  ev <- eigen(cor(X))$values
  nfact_kaiser <- sum(ev > 1)

  # ---------- 26.3 单因子性指标 ----------
  rmat <- cor(X)
  mean_r <- mean(rmat[upper.tri(rmat)])
  alpha20 <- tryCatch(psych::alpha(X, warnings = FALSE, check.keys = FALSE)$total$raw_alpha,
                      error = function(e) NA_real_)
  # 首因子解释比 = 第一特征值 / 总方差
  first_ratio <- ev[1] / sum(ev)

  # ---------- 26.2 EFA 载荷（1–6 因子，oblimin）----------
  maxload <- list()
  for (k in 1:6) {
    f <- tryCatch(
      psych::fa(X, nfactors = k, rotate = "oblimin", fm = "ml", warnings = FALSE),
      error = function(e) NULL
    )
    if (is.null(f)) next
    L <- as.data.frame(unclass(f$loadings))
    # 每题取绝对载荷最大的因子
    assign <- apply(abs(L), 1, which.max)
    maxv   <- apply(abs(L), 1, max)
    maxload[[as.character(k)]] <- data.frame(
      cohort = cohort_label, nfactors = k, item = rownames(L),
      max_loading = round(maxv, 3), assigned_factor = assign,
      row.names = NULL
    )
  }
  efa_tbl <- if (length(maxload)) do.call(rbind, maxload) else NULL

  # 载荷分布：各因子被分配给多少题（看是否形成 4 个清晰簇）
  cluster_tbl <- NULL
  if (!is.null(efa_tbl)) {
    cluster_tbl <- efa_tbl %>%
      group_by(cohort, nfactors, assigned_factor) %>%
      summarise(n_items = n(), mean_loading = round(mean(max_loading), 3), .groups = "drop") %>%
      arrange(nfactors, assigned_factor)
  }

  list(
    fit = data.frame(
      cohort = cohort_label, n = n,
      kaiser_factors = nfact_kaiser,
      parallel_factors = nfact_pa,
      eig1 = round(ev[1], 3), eig2 = round(ev[2], 3), eig3 = round(ev[3], 3),
      eig4 = round(ev[4], 3), eig5 = round(ev[5], 3),
      first_factor_ratio = round(first_ratio, 3),
      mean_inter_item_r = round(mean_r, 3),
      cronbach_20 = round(alpha20, 3),
      kmo = tryCatch(round(psych::KMO(X)$MSA, 3), error = function(e) NA_real_)
    ),
    efa = efa_tbl,
    cluster = cluster_tbl
  )
}

if (!all(ITEM_COLS %in% names(d))) {
  cat("  跳过：公开 L2 不含 hcf_5–hcf_24 原始题项；题目层 EFA 需要受控 L3。\n")
  .P3_MODULE_STATUS <- "skipped_controlled_l3_required"
  .P3_MODULE_REASON <- "需要受控 L3 的 hcf_5–hcf_24 原始题项"
} else {

res <- list()
res[["prolific"]] <- run_efa(d, "prolific", 462)
res <- res[!vapply(res, is.null, logical(1))]

if (!length(res)) {
  cat("  ⚠ 无可用队列，退出\n")
} else {
  # ---------- 写出 26.1：结构判据汇总 ----------
  fit_all <- do.call(rbind, lapply(res, `[[`, "fit"))
  write.csv(fit_all, file.path(OUT_DIR, "T26_item_structure.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n[写出] T26_item_structure.csv\n")
  print(fit_all)

  # ---------- 写出 26.2：EFA 载荷 ----------
  efa_all <- do.call(rbind, lapply(res, `[[`, "efa"))
  if (!is.null(efa_all)) {
    write.csv(efa_all, file.path(OUT_DIR, "T26_efa_loadings.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat("[写出] T26_efa_loadings.csv\n")
  }

  # ---------- 写出 26.3：因子-题目数分布（看有无清晰四簇）----------
  cl_all <- do.call(rbind, lapply(res, `[[`, "cluster"))
  if (!is.null(cl_all)) {
    write.csv(cl_all, file.path(OUT_DIR, "T26_efa_clusters.csv"),
              row.names = FALSE, fileEncoding = "UTF-8")
    cat("[写出] T26_efa_clusters.csv\n\n")
    cat("--- 各因子分配到的题数（nfactors=4 时应接近 5/5/5/5 才支持四因子）---\n")
    print(cl_all %>% filter(nfactors == 4))
  }

  # ---------- 结论口径 ----------
  cat("\n--- 结论口径 ---\n")
  for (nm in names(res)) {
    f <- res[[nm]]$fit
    cat(sprintf("%s：平行分析建议保留 %s 个因子；Kaiser(>1) = %s；首因子解释比 = %.3f；20 题 α = %.3f；平均题间 r = %.3f\n",
                nm, f$parallel_factors, f$kaiser_factors, f$first_factor_ratio,
                f$cronbach_20, f$mean_inter_item_r))
  }
  cat("\n  ⚠ 判读：若平行分析建议因子数 ≠ 4，或 nfactors=4 时题数分配严重不均，\n")
  cat("     则 20 题的**反身性四因子结构未获支持**——按 formative 立场如实报告，\n")
  cat("     并在论文中写明「以反身性模型检验形成性指数属模型误配」。\n")
}
}   # 原始题项齐全时执行
}   # ← 闭合 SKIP_JIANSHU 的 else 分支（2026-10-04：原为 quit() 提前退出）
