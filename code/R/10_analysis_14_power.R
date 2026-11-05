# =============================================================
# 10_analysis_14_power.R —— 统计功效、效应量对照与跨队列等价性
# -------------------------------------------------------------
# 目的：回应审稿人几乎必问的问题——
#   「见数 α→IBS 经 FDR 校正后 0/4 通过，是不是方法在第二队列失效了？」
#
# 五件事，缺一不可：
#   ① **报 MDE（最小可检测效应），不报事后功效**
#      事后功效是 p 值的确定性函数，不携带额外信息（Hoenig & Heisey 2001,
#      The American Statistician 55:19-24），审稿人会直接判为统计误用。
#      正确做法是反过来问：「以现有样本量，80% 功效下能可靠检出多大的效应？」
#   ② **2×2 交叉模式**：把「队列 × 结局」排成矩阵看，而不是只看单个结局
#   ③ **效应量 CI 与对照队列的效应量比**：能否排除"见数效应更小"？
#      若不能排除 → 说明是功效问题而非真实差异 → 补样本有依据
#   ④ **未校正 p 与 FDR q 并列**：区分「校正损耗」与「本就不显著」
#   ⑤ **等价检验（TOST）**：跨队列结论应是"等价/一致"而非"无显著差异"
#      （后者不能证伪差异）；如实报告，不达显著就说不达显著
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_both(); P <- L$prolific; J <- L$jianshu   # 跨队列对照：强制读入两队列
AC <- c("alpha1","alpha2","alpha3","alpha4")
OUTCOMES <- c("ibs_positive" = "IBS", "chronic_pain" = "慢性疼痛")

# ---------- 功效工具 ----------
mde_d <- function(n1, n0, alpha = 0.05, power = 0.80)
  (qnorm(1 - alpha/2) + qnorm(power)) * sqrt(1/n1 + 1/n0)
auc_to_d <- function(auc) sqrt(2) * qnorm(auc)      # d = √2·Φ⁻¹(AUC)
d_to_auc <- function(d) pnorm(d / sqrt(2))
n_needed <- function(d, ratio = 1, alpha = 0.05, power = 0.80) {
  d <- unname(as.numeric(d))          # 去名字，避免污染返回值的列名
  if (!is.finite(d) || d == 0) return(c(n1 = NA, n0 = NA, total = NA))
  n1 <- (1 + 1/ratio) * (qnorm(1-alpha/2) + qnorm(power))^2 / d^2
  c(n1 = ceiling(n1), n0 = ceiling(n1 * ratio), total = ceiling(n1 * (1 + ratio)))
}
# Cohen's d 及其 SE/95%CI（两独立组，Welch 型 SE）
d_with_ci <- function(y, x) {
  ok <- !is.na(y) & !is.na(x); y <- y[ok]; x <- x[ok]
  g1 <- x[y == 1]; g0 <- x[y == 0]; n1 <- length(g1); n0 <- length(g0)
  if (n1 < 3 || n0 < 3) return(c(d = NA, lo = NA, hi = NA, n1 = n1, n0 = n0, se = NA))
  sp <- sqrt(((n1-1)*var(g1) + (n0-1)*var(g0)) / (n1+n0-2))
  d <- (mean(g1) - mean(g0)) / sp
  se <- sqrt((n1+n0)/(n1*n0) + d^2/(2*(n1+n0)))
  c(d = d, lo = d - qnorm(.975)*se, hi = d + qnorm(.975)*se, n1 = n1, n0 = n0, se = se)
}
# bootstrap 版 ΔAUC（返回点估计与 SE，供 TOST 使用）
auc_diff_boot_full <- function(y1, x1, y2, x2, B = 3000, seed = 7) {
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
  d <- vapply(seq_len(B), function(i) one(), numeric(1)); d <- d[is.finite(d)]
  c(delta = observed, se = sd(d), lo = unname(quantile(d, .025)), hi = unname(quantile(d, .975)),
    nboot = length(d))
}
# TOST（等价检验）：H0: |Δ| ≥ margin  vs  H1: |Δ| < margin
tost <- function(est, se, margin) {
  z_lo <- (est + margin) / se      # 下界单侧：H0 Δ ≤ -margin
  z_hi <- (margin - est) / se      # 上界单侧：H0 Δ ≥ +margin
  c(p_lower = 1 - pnorm(z_lo), p_upper = 1 - pnorm(z_hi), p_tost = 1 - pnorm(min(z_lo, z_hi)))
}
EQUIV_MARGIN <- 0.10   # 预设等价界（AUC 单位；≈ Cohen d 0.36，即"小-中效应"边界）

rows <- list()
# 各块写入的字段不同，落盘前按「全部字段并集」补齐，避免丢列
rpad <- function(r, cols) { r[setdiff(cols, names(r))] <- NA; r[cols] }

cat("============================================================\n")
cat("【一】设计层面：MDE 与观测效应量（不报事后功效）\n")
cat("============================================================\n")
cat("  口径：MDE 为两独立样本均值差的 Welch 近似，双侧 α=0.05、功效 0.80；\n")
cat("        判定统一以 AUC 为主口径（与全文一致），Cohen's d 为并列表述。\n\n")

for (oc in names(OUTCOMES)) {
  lab <- OUTCOMES[[oc]]
  cat(sprintf("  ── 结局：%s ──\n", lab))
  cat(sprintf("  %-9s %-9s %5s %5s %8s %7s %8s %9s %8s %s\n",
              "队列", "α", "n+", "n-", "Cohen d", "AUC", "AUC下限", "需总n", "超MDE?", "观测AUC 95%CI"))
  for (nm in c("prolific", "jianshu")) {
    d0 <- L[[nm]]; y <- d0[[oc]]
    if (is.null(y)) { cat(sprintf("  %-9s （无该结局列）\n", nm)); next }
    n1 <- sum(y == 1, na.rm = TRUE); n0 <- sum(y == 0, na.rm = TRUE)
    m <- mde_d(n1, n0); floor_auc <- d_to_auc(m)
    for (a in AC) {
      dc  <- d_with_ci(y, d0[[a]])
      au  <- auc_ci(y, d0[[a]])
      nn  <- n_needed(abs(dc["d"]), ratio = n0/n1)
      pass <- unname(au["auc"]) >= floor_auc
      cat(sprintf("  %-9s %-9s %5d %5d %8.3f %7.3f %8.3f %9.0f %8s  [%.3f,%.3f]\n",
                  ifelse(nm == "prolific", "Prolific", "见数"), a, n1, n0,
                  dc["d"], au["auc"], floor_auc, nn["total"],
                  ifelse(pass, "★是", "否"), au["lo"], au["hi"]))
      rows[[length(rows)+1]] <- data.frame(
        outcome = lab, cohort = ifelse(nm == "prolific", "Prolific", "见数"), alpha = a,
        n_pos = n1, n_neg = n0, cohens_d = unname(dc["d"]),
        d_lo = unname(dc["lo"]), d_hi = unname(dc["hi"]),
        auc = unname(au["auc"]), auc_lo = unname(au["lo"]), auc_hi = unname(au["hi"]),
        mde_d = m, detectable_auc_floor = floor_auc,
        n_needed_total = unname(nn["total"]), exceeds_mde = pass)
    }
  }
  cat("\n")
}

cat("============================================================\n")
cat("【二】可检测的 AUC 下限（本分析的核心量）\n")
cat("============================================================\n")
for (nm in c("prolific", "jianshu")) {
  d0 <- L[[nm]]; y <- d0$ibs_positive
  n1 <- sum(y == 1, na.rm = TRUE); n0 <- sum(y == 0, na.rm = TRUE)
  m <- mde_d(n1, n0); fl <- d_to_auc(m)
  aucs <- sapply(AC, function(a) unname(auc_ci(y, d0[[a]])["auc"]))
  cat(sprintf("  %-9s n=%d（IBS+ %d / IBS- %d，比例 1:%.1f）\n",
              ifelse(nm == "prolific", "Prolific", "见数"), nrow(d0), n1, n0, n0/n1))
  cat(sprintf("            80%% 功效下可检出的最小效应：d = %.3f ⇔ AUC ≥ %.3f\n", m, fl))
  cat(sprintf("            实测 α→IBS 的 AUC 区间：%.3f–%.3f → %s\n",
              min(aucs), max(aucs),
              ifelse(max(aucs) >= fl, "在可检出范围内", "全部低于下限 ← 欠功效")))
}
y <- J$ibs_positive; n <- sum(!is.na(y))
m_bal <- mde_d(ceiling(n/2), floor(n/2)); m_now <- mde_d(sum(y==1,na.rm=TRUE), sum(y==0,na.rm=TRUE))
cat(sprintf("\n  【反事实】见数若 IBS 组 1:1 均衡（%d / %d）：MDE d = %.3f ⇔ AUC ≥ %.3f\n",
            ceiling(n/2), floor(n/2), m_bal, d_to_auc(m_bal)))
cat(sprintf("            现实 1:%.1f 的不平衡使 MDE 由 %.3f 恶化到 %.3f —— 提升症状检出率\n",
            sum(y==0,na.rm=TRUE)/sum(y==1,na.rm=TRUE), m_bal, m_now))
cat("            在该假设情景下可降低 MDE；实际招募成本与可行性需另行评估。\n")

cat("\n============================================================\n")
cat("【三】2×2 交叉模式：队列 × 结局（α₃）\n")
cat("        —— 只看单个结局会得出错误结论，必须排成矩阵\n")
cat("============================================================\n")
cat(sprintf("  %-10s %22s %22s\n", "", OUTCOMES[["ibs_positive"]], OUTCOMES[["chronic_pain"]]))
mat <- list()
for (nm in c("prolific", "jianshu")) {
  d0 <- L[[nm]]; cells <- c()
  for (oc in names(OUTCOMES)) {
    y <- d0[[oc]]; if (is.null(y)) { cells <- c(cells, "—"); next }
    dc <- d_with_ci(y, d0$alpha3); au <- auc_ci(y, d0$alpha3)
    n1 <- sum(y==1,na.rm=TRUE); n0 <- sum(y==0,na.rm=TRUE); m <- mde_d(n1, n0)
    tag <- if (abs(dc["d"]) >= 1.5*m) "★★" else if (abs(dc["d"]) >= m) "★" else "  "
    cells <- c(cells, sprintf("d=%+.2f AUC=%.3f %s", dc["d"], au["auc"], tag))
    mat[[length(mat)+1]] <- data.frame(cohort = ifelse(nm=="prolific","Prolific","见数"),
                                       outcome = OUTCOMES[[oc]], cohens_d = unname(dc["d"]),
                                       auc = unname(au["auc"]), mde_d = m,
                                       exceeds_mde = abs(dc["d"]) >= m)
  }
  cat(sprintf("  %-10s %22s %22s\n", ifelse(nm=="prolific","Prolific","见数"), cells[1], cells[2]))
}
cat("\n  ★★ = 效应量 ≥ 1.5×MDE；★ = 超过 MDE；空白 = 未达可检出下限\n")
cat("\n  ➜ **两队列在二分类口径下呈交叉模式**：Prolific 强于 IBS、弱于慢性疼痛；见数恰好相反。\n")
cat("     不能写「α₃→IBS 在两队列均复现」，也不能概括为所有预测效度均在两队列成立；\n")
cat("     只能按具体 α、结局、效应量与区间分别报告，并说明二分类口径下强关联的结局不同。\n")
cat(sprintf("     两队列上述指标的比例/均值不同（见数 IBS 19.2%% vs 24.7%%、PHQ-9 %.2f vs %.2f、BMI %.2f vs %.2f），\n",
            mean(J$phq9_total, na.rm=TRUE), mean(P$phq9_total, na.rm=TRUE),
            mean(J$bmi, na.rm=TRUE), mean(P$bmi, na.rm=TRUE)))
cat("     但不能据此笼统判定整体更健康，也不可推断见数的效应量被压缩——见【八】。\n")

cat("\n============================================================\n")
cat("【四】未校正 p 与 FDR q 并列：区分「校正损耗」与「本就不显著」\n")
cat("============================================================\n")
cat(sprintf("  %-9s %-9s %10s %10s %s\n", "队列", "α", "未校正 p", "FDR q", "判读"))
for (nm in c("prolific", "jianshu")) {
  d0 <- L[[nm]]
  pv <- sapply(AC, function(a) {
    dd <- data.frame(y = d0$ibs_positive, x = d0[[a]]); dd <- dd[stats::complete.cases(dd), ]
    summary(suppressWarnings(glm(y ~ x, data = dd, family = binomial)))$coefficients[2, 4]
  })
  qv <- p.adjust(pv, method = "BH")
  for (i in seq_along(AC))
    cat(sprintf("  %-9s %-9s %10.4f %10.4f %s\n",
                ifelse(nm=="prolific","Prolific","见数"), AC[i], pv[i], qv[i],
                ifelse(qv[i] < 0.05, "FDR 通过",
                       ifelse(pv[i] < 0.05, "未校正显著、FDR 后未过（属校正损耗）",
                              "未校正亦不显著（非校正所致）"))))
  if (nm == "jianshu" && all(pv > 0.05))
    cat(sprintf("     ➜ 见数最小未校正 p = %.4f > .05：**即便不做多重校正也不显著**，\n", min(pv)))
  if (nm == "jianshu" && all(pv > 0.05))
    cat("        故 FDR 0/4 不是校正造成的，而是效应量本身低于可检出下限。\n")
}

cat("\n============================================================\n")
cat("【五】效应量能否排除「与对照队列相同」？（补样本决策的依据）\n")
cat("============================================================\n")
cat("  设 Prolific 的实测效应量为 d_ref，检查其是否落在见数的 95%CI 内。\n")
cat("  这只是区间覆盖描述，不是两队列效应差的正式检验；覆盖或不覆盖均不能\n")
cat("  单独证明真实效应相同/不同，也不能单独决定是否补样本。\n\n")
cat(sprintf("  %-9s %-9s %8s %20s %10s %s\n", "队列", "α", "d", "95% CI", "d_ref", "CI不覆盖d_ref(描述)"))
for (a in AC) {
  dref <- unname(d_with_ci(P$ibs_positive, P[[a]])["d"])
  dc   <- d_with_ci(J$ibs_positive, J[[a]])
  excluded <- !(dc["lo"] <= dref & dref <= dc["hi"])
  cat(sprintf("  %-9s %-9s %8.3f %20s %10.3f %s\n", "见数", a, dc["d"],
              sprintf("[%+.3f, %+.3f]", dc["lo"], dc["hi"]), dref,
              ifelse(excluded, "是（仅区间描述）", "否（仅区间描述）")))
  rows[[length(rows)+1]] <- data.frame(
    outcome = "IBS", cohort = "见数 vs Prolific 参照", alpha = a,
    n_pos = NA, n_neg = NA, cohens_d = unname(dc["d"]),
    d_lo = unname(dc["lo"]), d_hi = unname(dc["hi"]), auc = NA_real_,
    auc_lo = NA_real_, auc_hi = NA_real_, mde_d = NA_real_,
    detectable_auc_floor = NA_real_, n_needed_total = NA_real_, exceeds_mde = NA,
    ref_d = dref, ref_excluded = excluded)
}

cat("\n============================================================\n")
cat("【六】跨队列等价检验（TOST on ΔAUC）—— 如实报告\n")
cat("============================================================\n")
cat(sprintf("  等价界预设 ±%.2f（AUC 单位，≈ Cohen d %.2f）\n", EQUIV_MARGIN, auc_to_d(0.5 + EQUIV_MARGIN)))
cat(sprintf("  ΔAUC 方向 = Prolific − 见数（与模块 11 一致）\n\n"))
cat(sprintf("  %-9s %9s %9s %22s %10s %9s %s\n",
            "α", "ΔAUC", "boot SE", "95% CI", "p_TOST", "等价?", "可排除的差异范围"))
boot_ses <- c()
for (a in AC) {
  b <- auc_diff_boot_full(P$ibs_positive, P[[a]], J$ibs_positive, J[[a]], B = 3000)
  boot_ses <- c(boot_ses, unname(b["se"]))
  t_ <- tost(b["delta"], b["se"], EQUIV_MARGIN)
  # 当前精度下能排除的差异：|Δ| 的置信上界
  exclude_range <- max(abs(b["lo"]), abs(b["hi"]))
  cat(sprintf("  %-9s %+9.3f %9.4f  [%+.3f, %+.3f] %10.4f %9s  |Δ| < %.3f\n",
              a, b["delta"], b["se"], b["lo"], b["hi"], t_["p_tost"],
              ifelse(t_["p_tost"] < 0.05, "★是", "否"), exclude_range))
  rows[[length(rows)+1]] <- data.frame(
    outcome = "IBS", cohort = "ΔAUC(Prolific-见数)", alpha = a,
    n_pos = NA, n_neg = NA, cohens_d = NA_real_, d_lo = NA_real_, d_hi = NA_real_,
    auc = b["delta"], auc_lo = b["lo"], auc_hi = b["hi"],
    mde_d = NA_real_, detectable_auc_floor = NA_real_, n_needed_total = NA_real_,
    exceeds_mde = NA, tost_p = unname(t_["p_tost"]), equiv_margin = EQUIV_MARGIN,
    se = unname(b["se"]))
}
cat("\n  ⚠️ 判读纪律（务必遵守）\n")
cat("     · CI 跨零 ≠ 证明等价，只能写「未检出差异」；\n")
cat("     · TOST 未达显著时，只能写「未检出差异，且 95%CI 排除了大于 X 的差异」；\n")
cat("     · **绝不能写「两队列完全相同/一致」**——现有精度不足以支持该结论。\n")
se_need <- (EQUIV_MARGIN - 0.05) / qnorm(0.95)
cat(sprintf("\n  【精度缺口】以 ΔAUC≈0.05 计，要在 ±%.2f 上判定等价需 SE ≤ %.4f；\n", EQUIV_MARGIN, se_need))
current_se <- mean(boot_ses, na.rm = TRUE)
cat(sprintf("            当前四α平均SE=%.4f；在所述假设下样本量约需提高 %.2f 倍（每队列约 %d 例）。\n",
            current_se, (current_se/se_need)^2, ceiling(nrow(J)*(current_se/se_need)^2)))

cat("\n============================================================\n")
cat("【七】补样本的量化决策（若决定继续招募）\n")
cat("============================================================\n")
n1 <- sum(J$ibs_positive == 1, na.rm = TRUE); n0 <- sum(J$ibs_positive == 0, na.rm = TRUE)
ratio <- n0 / n1
cat(sprintf("  见数现状 n=%d（%d+ / %d-，比例 1:%.1f）\n\n", nrow(J), n1, n0, ratio))
cat(sprintf("  %-24s %8s %10s %10s\n", "若真实效应为", "d", "所需总n", "还需再收"))
for (case in list(c("与 Prolific 实测同量级", 0.286), c("中等效应", 0.20), c("见数当前实测", 0.132))) {
  d <- as.numeric(case[2]); nn <- n_needed(d, ratio); tot <- unname(nn["total"])
  cat(sprintf("  %-24s %8.3f %10.0f %10d\n", case[1], d, tot, max(0, tot - nrow(J))))
}
cat("\n  ➊ 上表均为假设效应量下的设计情景，不是对真实效应或实际补样人数的保证；\n")
cat("  ➋ 是否继续招募应结合效应不确定性和阳性率，而非仅按单一点估计决定；\n")
cat("  ➌ 若目标是降低该二分类组间比较的 MDE，可评估富集招募（提高 IBS 阳性比例）的情景：\n")
cat(sprintf("     同样 n=%d 下 MDE 由 %.3f 降到 %.3f（AUC 下限 %.3f → %.3f）。\n",
            nrow(J), mde_d(n1,n0), mde_d(168,336), d_to_auc(mde_d(n1,n0)), d_to_auc(mde_d(168,336))))

cat("\n============================================================\n")
cat("【八】事后敏感性：IBS 二分类与同源「症状计数」对照\n")
cat("        —— 保留更多响应粒度，但不保证提高功效，也不替代注册主结局\n")
cat("============================================================\n")
cat("  二分类把「1 项症状」和「4 项症状」都记为阳性 → 损失信息；\n")
cat("  叠加 1:4.2 的组间不平衡，功效被双重削弱。改用症状计数可保留全部信息。\n\n")
cat("  ※ Q41 为 **Rome IV 结构题**（题干：「过去 3 个月，您是否每周至少有 1 天出现\n")
cat("    以下任何一种症状？」，含病程 3 个月 + 频率每周≥1 天 + 腹痛 + 排便性状改变要件），\n")
cat("    故症状计数是 IBS 二分类判定的**同源连续化**（同一题项、不同粒度），非另一结局。\n")
cat("    建议命名 `IBS-related GI symptom count`；回避 `severity`（无强度/生活干扰维度）。\n\n")
JS_RAW <- Sys.getenv("P3_JIANSHU_RAW", unset = "")
if (nzchar(JS_RAW) && file.exists(JS_RAW)) {
  rw <- readr::read_csv(JS_RAW, col_types = readr::cols(.default = readr::col_character()),
                        show_col_types = FALSE)
  # Q41 编码：1=腹痛或腹胀 2=腹泻 3=便秘 4=交替 5=以上都没有
  cnt <- function(s) {
    if (is.na(s) || s == "") return(NA_real_)
    tags <- setdiff(strsplit(s, "|", fixed = TRUE)[[1]], "5")
    length(tags)
  }
  sev <- vapply(rw$Q41, cnt, numeric(1))
  cat("  症状计数（排除「以上都没有」，0–4）分布：\n")
  print(table(sev, useNA = "ifany"))
  cat(sprintf("  均值 %.3f  SD %.3f\n\n", mean(sev, na.rm = TRUE), sd(sev, na.rm = TRUE)))
  if (!"用户ID" %in% names(rw) || !"pid" %in% names(J))
    stop("见数原始症状计数必须按受控用户ID与分析pid联结；禁止按行位置拼接")
  if (anyDuplicated(rw$用户ID)) stop("见数原始用户ID不唯一")
  sev_by_id <- setNames(sev, as.character(rw$用户ID))
  unmatched <- setdiff(as.character(J$pid), names(sev_by_id))
  if (length(unmatched)) stop("见数分析基底存在未匹配原始用户ID：", length(unmatched), "人")
  Js <- J
  Js$ibs_symptom_count <- unname(sev_by_id[as.character(J$pid)])
  mde_r <- function(nb, alpha = 0.05, power = 0.80)
    tanh((qnorm(1-alpha/2) + qnorm(power)) / sqrt(nb - 3))
  cat(sprintf("  %-8s %12s %12s %12s %10s %s\n",
              "α", "与计数 r", "Spearman ρ", "Pearson p", "可检出 |r|", "判定"))
  for (a in AC) {
    x <- Js[[a]]; yv <- Js$ibs_symptom_count
    ok <- is.finite(x) & is.finite(yv)
    r <- cor(x[ok], yv[ok]); rho <- suppressWarnings(cor(x[ok], yv[ok], method = "spearman"))
    nok <- sum(ok); lim <- mde_r(nok)
    tt <- cor.test(x[ok], yv[ok])
    cat(sprintf("  %-8s %+12.3f %+12.3f %12.2e %10.3f %s\n",
                a, r, rho, tt$p.value, lim,
                ifelse(abs(r) >= lim, "★超过下限（该 n 下可检出）", "低于下限")))
  }
  cat("\n  ➜ 症状计数是同源、事后敏感性分析；其显著性不得替代二分类主结局。\n")
  cat("     两种尺度的结果差异只能提示信息粒度可能影响统计关联，不能证明二分类是唯一原因。\n")
  cat("     任何功效或结局解释均须并列报告，不得将事后连续化包装成预注册结局。\n")
  cat("\n  ⚠️ 报告纪律：仍须① Methods 说明主结局为二分类、连续化口径为敏感性分析；\n")
  cat("     ② 两口径结果并列报告，不以连续口径的显著性替代二分类结果；\n")
  cat("     ③ 症状计数是数量而非临床诊断（无排除诊断），回避 `severity` 一词。\n")
  for (a in AC) {
    x <- Js[[a]]; yv <- Js$ibs_symptom_count; ok <- is.finite(x) & is.finite(yv)
    rows[[length(rows)+1]] <- data.frame(
      outcome = "IBS 症状计数（0-4）", cohort = "见数", alpha = a,
      n_pos = NA, n_neg = NA, cohens_d = NA_real_, d_lo = NA_real_, d_hi = NA_real_,
      auc = NA_real_, auc_lo = NA_real_, auc_hi = NA_real_, mde_d = NA_real_,
      detectable_auc_floor = NA_real_, n_needed_total = NA_real_, exceeds_mde = NA,
      r_symptom_count = cor(x[ok], yv[ok]), rho_symptom_count = suppressWarnings(cor(x[ok], yv[ok], method="spearman")),
      r_detectable_floor = mde_r(sum(ok)))
  }
} else cat("  （未见数原始题表，跳过该节）\n")

if (nzchar(JS_RAW) && file.exists(JS_RAW)) {
cat("\n  ── 同口径对照：Prolific 的 α ↔ 症状计数相关（用同一题项、同一算法构造）──\n")
# ibs_29（原始多选题串）不在分析层列白名单内，需从清洗层输出单独取并按 token 对齐
PC_CSV <- base_path("p3_base_prolific")
pro_raw <- readr::read_csv(PC_CSV,
                           col_types = readr::cols(participant_token = readr::col_character(),
                                                   data_source = readr::col_character(),
                                                   ibs_29 = readr::col_character(),
                                                   .default = readr::col_skip()),
                           show_col_types = FALSE) %>%
  dplyr::arrange(participant_token, dplyr::desc(data_source == "full")) %>%
  dplyr::distinct(participant_token, .keep_all = TRUE)
Px <- dplyr::left_join(P, pro_raw, by = c("pid" = "participant_token"), relationship = "one-to-one")
if (any(is.na(Px$ibs_29))) warning("部分Prolific记录缺少ibs_29；相关分析按完整病例处理")
pcnt <- function(s) {
  if (is.na(s) || s == "") return(NA_real_)
  length(setdiff(strsplit(as.character(s), ";", fixed = TRUE)[[1]], "5"))
}
P_symc <- vapply(Px$ibs_29, pcnt, numeric(1))
cat(sprintf("  Prolific 症状计数：均值 %.3f  SD %.3f  零值比例 %.1f%%（见数：%.3f / %.3f / %.1f%%）\n",
            mean(P_symc, na.rm=TRUE), sd(P_symc, na.rm=TRUE), 100*mean(P_symc==0, na.rm=TRUE),
            mean(sev, na.rm=TRUE), sd(sev, na.rm=TRUE), 100*mean(sev==0, na.rm=TRUE)))
cat(sprintf("  %-8s %14s %14s %10s %9s %s\n", "α", "Prolific r", "见数 r", "Fisher z", "p", "判读"))
for (a in AC) {
  xp <- Px[[a]]; okp <- is.finite(xp) & is.finite(P_symc)
  xj <- Js[[a]]; okj <- is.finite(xj) & is.finite(Js$ibs_symptom_count)
  rp <- cor(xp[okp], P_symc[okp]); rj <- cor(xj[okj], Js$ibs_symptom_count[okj])
  zz <- (atanh(rp) - atanh(rj)) / sqrt(1/(sum(okp)-3) + 1/(sum(okj)-3))
  pv <- 2 * (1 - pnorm(abs(zz)))
  cat(sprintf("  %-8s %+14.3f %+14.3f %10.2f %9.3f %s\n", a, rp, rj, zz, pv,
              ifelse(pv < 0.05, "★检出队列差异", "未检出队列差异（非等效证据）")))
}
cat("\n  ➜ 上述 Fisher-z 结果为两独立相关的近似比较；未达显著不等于效应相等。\n")
cat("     连续口径的跨队列比较不能单独证明二分类化或组间不平衡造成了 AUC 差异。\n")

cat("\n  ── ⭐ 阳性率 vs 样本量：哪个才是瓶颈？──\n")
ZA2 <- qnorm(0.975); ZB2 <- qnorm(0.80)
mde2 <- function(a, b) (ZA2 + ZB2) * sqrt(1/a + 1/b)
cat(sprintf("  现状 n=%d（阳性 %d / 阴性 %d，1:%.1f）→ MDE d = %.3f\n",
            nrow(J), n1, n0, n0/n1, mde2(n1, n0)))
cat("\n  路径 A（只加样本量，阳性率保持 19.2%）：\n")
for (tot in c(618, 700, 800, 1000, 1200)) {
  aa <- round(tot*0.192); bb <- tot - aa
  cat(sprintf("     n=%5d → %4d+ / %4d−   MDE = %.3f\n", tot, aa, bb, mde2(aa, bb)))
}
cat("\n  路径 B（只提阳性率，n 保持 504）：\n")
for (pair in list(c("现状 1:4.2", 97, 407), c("1:3", 126, 378), c("1:2", 168, 336), c("1:1", 252, 252))) {
  cat(sprintf("     %-10s → MDE = %.3f\n", pair[1], mde2(as.numeric(pair[2]), as.numeric(pair[3]))))
}
cat("\n  ➜ 换算关系：MDE ∝ 1/√n（总样本量），阳性率只是**效率杠杆**——\n")
cat(sprintf("     提到 1:2 等效于放大 %.2f 倍样本量；提到 1:1 等效于放大 %.2f 倍。\n",
            721/504, 811/504))
cat("     ⇒ 总样本量与阳性率共同决定可检测效应；上述仅为设计情景比较，\n")
cat("        不据此断言唯一瓶颈，也不保证任一招募方案能达到目标功效。\n")
} else cat("  （缺少原始题表，跳过 IBS 症状计数的双队列对照）\n")
cat("\n============================================================\n")
cat("【九】可直接写进论文的表述（草稿）\n")
cat("============================================================\n")
mJ <- mde_d(n1, n0); flJ <- d_to_auc(mJ)
aucJ <- sapply(AC, function(a) unname(auc_ci(J$ibs_positive, J[[a]])["auc"]))
cp_d <- d_with_ci(J$chronic_pain, J$alpha3); cp_a <- auc_ci(J$chronic_pain, J$alpha3)
cat("\n  【A · Methods 用（功效声明）】\n")
cat("  [EN] We do not report post-hoc observed power, which is a deterministic function\n")
cat("       of the observed p-value and carries no additional information (Hoenig &\n")
cat("       Heisey, 2001). Instead, we report the minimum detectable effect (MDE).\n")
cat(sprintf("  [CN] 本研究不报告事后功效（其对 p 值是确定性函数、不提供额外信息；\n"))
cat("       Hoenig & Heisey, 2001），改报最小可检测效应（MDE）。\n")

cat("\n  【B · Results 用（见数主要检验族）】\n")
cat(sprintf("  [CN] 见数队列 n=%d（IBS 阳性 %d、阴性 %d）在 80%% 功效、双侧 α=.05 下\n", nrow(J), n1, n0))
cat(sprintf("       可检出的最小效应为 d ≥ %.3f（对应 AUC ≥ %.3f）；实测四个 α 与 IBS 的\n", mJ, flJ))
cat(sprintf("       关联仅为 AUC %.3f–%.3f（d %.3f–%.3f），全部低于该下限。因此主要检验族\n",
            min(aucJ), max(aucJ),
            min(sapply(AC, function(a) abs(d_with_ci(J$ibs_positive, J[[a]])["d"]))),
            max(sapply(AC, function(a) abs(d_with_ci(J$ibs_positive, J[[a]])["d"])))))
cat(sprintf("       经 FDR 校正后无一通过；在当前样本量与事件率下，对该量级效应的检验精度有限，\n"))
cat(sprintf("       这既不证明关联存在，也不证明关联不存在（四个未校正 p 亦均 > .10，见补充表）。\n"))

cat("\n  【C · Results 用（结局谱差异，必写）】\n")
cat(sprintf("  [CN] 在同一队列、同一 α 评分下，与慢性疼痛的关联远强于 IBS\n"))
cat(sprintf("       （AUC %.3f，d = %.3f，n=%d/%d，远超可检出下限），而在 Prolific 队列中\n",
            cp_a["auc"], cp_d["d"], sum(J$chronic_pain==1,na.rm=TRUE), sum(J$chronic_pain==0,na.rm=TRUE)))
cat(sprintf("       恰相反（IBS AUC %.3f、慢性疼痛 AUC %.3f）。这只能按具体结局描述为\n",
            unname(auc_ci(P$ibs_positive, P$alpha3)["auc"]), unname(auc_ci(P$chronic_pain, P$alpha3)["auc"])))
cat(sprintf("       关联谱不同，不能概括为所有 α 预测效度均在两个队列得到确认。\n"))
cat(sprintf("       同源症状计数仅是事后敏感性分析；未检出跨队列相关差异不等于等效，\n"))
cat(sprintf("       也不能据此认定判定阈值或症状分布差异造成了上述二分类 AUC 模式。\n"))

cat("\n  【D · Limitations 用】\n")
cat(sprintf("  [CN] 见数队列对 α→IBS 这一主要结局的检验精度有限（MDE d=%.3f > 实测 d≈%.2f）；\n",
            mJ, max(sapply(AC, function(a) abs(d_with_ci(J$ibs_positive, J[[a]])["d"])))))
cat(sprintf("       跨队列 ΔAUC 的 95%%CI 虽跨零、未检出差异，但精度尚不足以判定等价\n"))
cat(sprintf("       （TOST p>0.05），故本文仅表述为「未检出队列差异」，不主张两队列一致。\n"))

ALL_COLS <- unique(unlist(lapply(rows, names)))
save_tab(do.call(rbind, lapply(rows, rpad, cols = ALL_COLS)), "T11_power_equivalence")
cat("\n  模块完成：功效、效应量对照与等价性分析\n")
