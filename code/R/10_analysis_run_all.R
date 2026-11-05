# =============================================================
# 10_analysis_run_all.R —— 分析层总驱动
# 依次执行 01–19 号模块，全部控制台输出落盘到 <结果目录>/analysis_log.txt
# 用法： Rscript 10_analysis_run_all.R
# 结果目录由 P3_COHORT 决定：
#   prolific → <项目>/数据分析/分析结果
#   jianshu  → <项目>/见数/分析结果
#   both     → <项目>/output
# 可用 P3_OUT_DIR 显式指定隔离输出目录；总入口清单/日志与各模块结果表共用该目录
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
PROJ <- Sys.getenv("P3_ROOT", unset = getwd())
BASE <- Sys.getenv("P3_CODE_DIR", unset = file.path(PROJ, "code", "R"))
P3_COHORT <- Sys.getenv("P3_COHORT", "both")     # 公开总入口默认双队列；单队列运行需显式设置环境变量
OUT <- Sys.getenv("P3_OUT_DIR", unset = switch(tolower(P3_COHORT),
               prolific = file.path(PROJ, "数据分析", "分析结果"),
               jianshu  = file.path(PROJ, "见数", "分析结果"),
               file.path(PROJ, "output")))
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUT, "analysis_log.txt")
con  <- file(logf, open = "wt", encoding = "UTF-8")
sink(con, split = TRUE)
cat("运行模式 P3_COHORT =", P3_COHORT, "｜ 结果输出目录：", OUT, "\n")

scripts <- c("10_analysis_01_descriptives.R","10_analysis_02_alpha_structure.R",
             "10_analysis_03_causal_chain.R","10_analysis_04_predictive.R",
             "10_analysis_05_ibs_dose_response.R","10_analysis_06_clinical_utility.R",
             "10_analysis_07_hidden_distress.R","10_analysis_09_carrier.R",
             "10_analysis_11_cross_cohort.R","10_analysis_12_convergent_semantic.R",
             "10_analysis_14_power.R",
             "10_analysis_15_screening_roc.R","10_analysis_16_hierarchical_regression.R",
             "10_analysis_17_screening_equivalence.R",             "10_analysis_18_anchor_vs_algorithm.R",
             # ---- 20–24：零成本证据补强（2026-09-30 新增）----
             "10_analysis_20_external_validation.R",     # 跨队列切点外部验证
             "10_analysis_21_acceptability.R",           # 可接受性/可行性/时间成本
             "10_analysis_22_incremental_mtmm.R",        # 增量效度 + MTMM（仅 Prolific）
             "10_analysis_23_clinical_operability.R",    # PPV/NPV 随患病率 + 已知组效度
             "10_analysis_24_item_equivalence.R",        # 跨队列题目等价性
             # ---- 25：复现核定与扩展 ----
             "10_analysis_25_replication_extension.R",   # 肥胖OR/BMI/数值锚点信度/LPA/外部分布
             # ---- 26：题目结构 EFA（2026-10-03 新增）----
             "10_analysis_26_item_structure_efa.R",      # 20题因子结构（无预设归属，不触碰映射）
             # ---- 27：H8 自杀意念条目（2026-10-04 新增）----
             "10_analysis_27_suicidal_ideation.R",     # α₁×PHQ-9 第9条（双队列）
             # ---- 32：设计层样本量 ----
             "10_analysis_32_design_power.R",             # E3a／E3b 最小样本量论证（设计层）
             "10_analysis_34_disease_or.R",              # 七疾病 OR（补齐 A2/C11/C15：IBS 外五病 α~慢病）
             # ---- 35：跨队列结局对照（2026-10-16 新增，补齐用户要求的肥胖/慢性痛/IBS剂量反应/NNT/切点/U型/决策树 跨队列比）----
             "10_analysis_35_cross_cohort_outcomes.R",   # load_both() 强制双队列，统一重算跨队列对比
             # ---- 36：7 因子主分析（2026-10-26 新增）----
             #   7 因子 = α₁/α₂/α₄ + 载体 a/b/c/d；统一使用单一固定管线读数
             #   （直接采用后台 v1.1 参数导出的 scores_* 基线值）
             "10_analysis_36_seven_factors.R",
             # ---- 43：E3b 注册 C1/C2/C3（公开输出层）----
             "10_analysis_43_e3b_registered.R")
failed <- c()
skipped <- c()
partial <- c()
partial_controlled <- c()
public_header <- names(read.csv(file.path(PROJ, "code", "data", "derived", "p3_base_prolific.csv"),
                               nrows = 1, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
public_l2_mode <- !any(c("pid", "hcf_5", "phq_36", "anchor1") %in% public_header)
partial_reasons <- c(
  "10_analysis_21_acceptability.R" = "公开 L2 不含实测时长/注意力字段；仅公开分布与缺失描述可运行",
  "10_analysis_24_item_equivalence.R" = "公开 L2 不含题目级输入；仅 α/H 分布比较可运行",
  "10_analysis_25_replication_extension.R" = "公开 L2 不含锚点题；锚定信度未运行，其余扩展分析已完成"
)
module_reasons <- setNames(rep("", length(scripts)), scripts)
# 记录总入口顶层变量，供运行状态核验使用。
runall_own <- ls(.GlobalEnv, all.names = TRUE)
# 保存平台目录副本，供防污染守卫比对（模块若覆盖 BASE / OUT 会被自动恢复）
.PLAT_BASE <- BASE
.PLAT_OUT  <- OUT   # 多个模块（07/08/10/13/20）把 OUT 当作本地 accumulator 复用，会污染输出目录
for (.module_script in scripts) {
  cat("\n\n############################################################\n")
  cat("####  ", .module_script, "\n")
  cat("############################################################\n")
  .P3_MODULE_STATUS <- "success"
  .P3_MODULE_REASON <- ""
  ok <- tryCatch({
    # ⚠️ 逐模块重设 locale：部分模块用 readr 直接读中文路径 CSV，
    #    批处理语境下 locale 非 UTF-8 会报「无法打开链结」。
    suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
    source(file.path(BASE, .module_script)); TRUE
  },
                 error = function(e) {
                   cat("!!! 该脚本运行失败 !!!\n错误信息：", conditionMessage(e), "\n")
                   .P3_MODULE_STATUS <<- "failed"
                   .P3_MODULE_REASON <<- conditionMessage(e)
                   FALSE
                 })
  # ⚠️ 防污染守卫：模块若在顶层给 BASE 赋值，会覆盖 run_all 的平台目录，
  #    导致其后所有模块 source 到错误路径（症状：14 号起连续报「无法打开链结」）。
  #    2026-10-04 已修 13 号；此处再加校验，一旦被改立即恢复并告警。
  if (!identical(BASE, .PLAT_BASE)) {
    cat(sprintf("⚠️ 警告：模块 %s 污染了全局 BASE，已自动恢复\n", .module_script))
    cat(sprintf("   被改为: %s\n   应为  : %s\n", BASE, .PLAT_BASE))
    assign("BASE", .PLAT_BASE, envir = .GlobalEnv)
  }
  # OUT 污染守卫：模块若复用 OUT 作累加器，恢复总入口输出目录。
  if (!identical(OUT, .PLAT_OUT)) {
    cat(sprintf("⚠️ 警告：模块 %s 污染了全局 OUT（输出目录），已自动恢复\n", .module_script))
    assign("OUT", .PLAT_OUT, envir = .GlobalEnv)
  }
  if (!ok || identical(.P3_MODULE_STATUS, "failed")) {
    failed <- c(failed, .module_script)
    module_reasons[.module_script] <- .P3_MODULE_REASON
  } else if (identical(.P3_MODULE_STATUS, "skipped_controlled_l3_required")) {
    skipped <- c(skipped, .module_script)
    module_reasons[.module_script] <- .P3_MODULE_REASON
    cat(">>> 受控 L3 专属模块显式跳过：", .P3_MODULE_REASON, "\n")
  } else if (public_l2_mode && .module_script %in% names(partial_reasons)) {
    partial <- c(partial, .module_script)
    module_reasons[.module_script] <- partial_reasons[[.module_script]]
    cat(">>> 公开 L2 部分完成：", module_reasons[[.module_script]], "\n")
  } else if (identical(.P3_MODULE_STATUS, "partial_controlled_l3")) {
    partial_controlled <- c(partial_controlled, .module_script)
    module_reasons[.module_script] <- .P3_MODULE_REASON
    cat(">>> 受控 L3 部分完成：", .P3_MODULE_REASON, "\n")
  }
}

module_status <- ifelse(scripts %in% failed, "failed",
                 ifelse(scripts %in% skipped, "skipped_controlled_l3_required",
                 ifelse(scripts %in% partial, "partial_public_l2",
                 ifelse(scripts %in% partial_controlled, "partial_controlled_l3", "success"))))
manifest <- data.frame(module = scripts, status = module_status,
                       reason = unname(module_reasons[scripts]), stringsAsFactors = FALSE)
write.csv(manifest, file.path(OUT, "run_all_manifest.csv"),
          row.names = FALSE, fileEncoding = "UTF-8", na = "")
cat("运行清单：", file.path(OUT, "run_all_manifest.csv"), "\n")
sink(); close(con)
cat("全部完成，日志：", logf, "\n")
if (length(failed)) {
  cat("⚠ 以下脚本失败：", paste(failed, collapse = ", "), "\n")
  quit(save = "no", status = 1)
} else {
  cat(sprintf("✓ %d 个模块完整成功；%d 个公开 L2 部分完成；%d 个受控 L3 部分完成；%d 个受控 L3 专属模块显式跳过；0 个失败\n",
              length(scripts) - length(skipped) - length(partial) - length(partial_controlled),
              length(partial), length(partial_controlled), length(skipped)))
}
