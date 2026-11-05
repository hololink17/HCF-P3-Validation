# =============================================================
# 10_analysis_07_hidden_distress.R —— R1/R5 分型 + 潜在痛苦（论文头号卖点）
# 定义（与 Prolific 原始 代码1 逐字一致，见 R/00_shared_rules.R::add_alpha_typing）：
#   分型  = 四个 α 各按本队列 P33/P67 分低/中/高，四者同档才归类（否则「其他」）
#   潜在痛苦 = PHQ-9 正常（<10）但 α 归类为「耗竭型」→ 传统量表"测不到"的人群
# 对应已发表论文：
#   · Condition 2（现象学状态定义）—— 四类状态独立于测量工具
#   · Condition 4 / 增量效度（回应"传统量表原理上测不到"）
# 两队列同规则；见数无「其他」以外的后端分型，一律以本规则重建
# =============================================================
# 🔴 2026-10-06 口径修正与扩展（旧稿数值一律作废，不在此处复写）：
#   ① 7.5 节原**只 cat 到日志、未落任何表** → 现落 T7_incremental_vs_phq9.csv，
#      并补 bootstrap 95% CI（set.seed 固定 → 可复现）。
#   ② ⚠️ 对照口径不可混写：本节对照 ＝ **单一症状量表（PHQ-9）**，样本为全队列；
#      模块 22 的 T22_incremental_validity 对照 ＝ **四金标准组合**，样本为
#      **Prolific 完整版 n=142 且仅单队列**。两者是不同对照、不同样本，
#      引用时必须同句并列（"对照成熟量表组合→未达显著；对照单一症状量表→显著且双队列复现"），
#      **禁止合并成一句、禁止互相比大小**。
#   ③ 各结局**各自 drop_na**。旧稿两个结局共用一个 drop_na 列表，
#      会让二值结局（IBS）的缺失反过来削减连续结局（GAD-7）的样本。
#   ④ 新增 7.6 节「双正常亚组」（PHQ-9<10 且 GAD-7<10）：把"阈下识别"推进到
#      两个症状量表均低于临床阈值的更严格子集。落 T7_doubly_normal_split.csv。
#   ⑤ 🔴 7.6 实测**证伪**了"结构分层与症状解耦"这一预设：组间 PHQ-9/GAD-7 差
#      未压到 0（Prolific d≈1.0，见数 d≈1.6–1.9）。故本表只支持"阈下亚临床梯度"，
#      **不得**写成"测到完全无症状者的隐藏问题"。
#   ⑥ 7.3 节的组间对比此前同样**只 print 未落表**（E7 缺可追溯表源）→
#      现落 T7_potential_distress_vs_healthy.csv，并补 d 的 bootstrap 95% CI。
#   ⑦ H 值与分型同源，其组间差仅作描述，不作为独立外部效度证据。
#      一类），两表均加"［分组定义变量·循环，仅供核对］"标记，不得作为结果引用。
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts()
OUT <- list(); GT <- list(); PD <- list(); DN <- list(); IV <- list()

# ---------- bootstrap 工具（种子固定 → 结果可复现）----------
BS_B <- 2000; BS_SEED <- 20261006

# 独立两组 Cohen's d 的 bootstrap 95% CI（各组内重抽样）
cohen_d_ci <- function(x, y, B = BS_B, seed = BS_SEED) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 3 || length(y) < 3) return(c(NA_real_, NA_real_))
  set.seed(seed)
  v <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    xi <- sample(x, length(x), replace = TRUE); yi <- sample(y, length(y), replace = TRUE)
    s <- sqrt((var(xi) + var(yi)) / 2)
    if (is.finite(s) && s > 0) v[b] <- (mean(xi) - mean(yi)) / s
  }
  quantile(v, c(.025, .975), na.rm = TRUE)
}

# α₃ 增量（在 phq9_total 之上）的 bootstrap 95% CI：连续结局取 ΔR²，二值结局取 ΔAUC
bs_delta <- function(dat, y, kind, B = BS_B, seed = BS_SEED) {
  set.seed(seed)
  v <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    i <- sample.int(nrow(dat), nrow(dat), replace = TRUE)
    dd <- dat[i, , drop = FALSE]
    if (kind == "cont") {
      m0 <- lm(reformulate("phq9_total", y), data = dd)
      m1 <- lm(reformulate(c("phq9_total", "alpha3"), y), data = dd)
      v[b] <- summary(m1)$r.squared - summary(m0)$r.squared
    } else {
      if (length(unique(stats::na.omit(dd[[y]]))) < 2) next
      m0 <- suppressWarnings(glm(reformulate("phq9_total", y), data = dd, family = binomial))
      m1 <- suppressWarnings(glm(reformulate(c("phq9_total", "alpha3"), y), data = dd, family = binomial))
      v[b] <- as.numeric(pROC::auc(pROC::roc(dd[[y]], predict(m1, type = "response"), quiet = TRUE))) -
              as.numeric(pROC::auc(pROC::roc(dd[[y]], predict(m0, type = "response"), quiet = TRUE)))
    }
  }
  quantile(v, c(.025, .975), na.rm = TRUE)
}

for (nm in names(L)) {
  d <- L[[nm]]
  cat(sprintf("\n=== %s：分型与潜在痛苦（N=%d）===\n", nm, nrow(d)))

  # ---- 7.1 四类分型（含「其他」）----
  g <- d %>% group_by(alpha_group) %>%
    summarise(n = n(), H = mean(H_score, na.rm = TRUE), phq9 = mean(phq9_total, na.rm = TRUE),
              gad7 = mean(gad7_total, na.rm = TRUE),
              ibs = mean(ibs_positive, na.rm = TRUE) * 100,
              dep_rate = mean(phq9_total >= 10, na.rm = TRUE) * 100, .groups = "drop") %>%
    arrange(factor(alpha_group, levels = c("耗竭型","均衡型","高意识型","其他")))
  cat("  分型（profile 型，含「其他」）:\n"); print(as.data.frame(g %>% mutate(H = round(H,3), phq9 = round(phq9,1), gad7 = round(gad7,1), ibs = round(ibs,1), dep_rate = round(dep_rate,1))))

  # ---- 7.2 潜在痛苦 / 真正健康 ----
  norm <- d %>% filter(phq9_total < 10)
  n_pd <- sum(norm$potential_depression); n_hc <- sum(norm$healthy_control)
  cat(sprintf("  PHQ-9 正常者 n=%d；潜在痛苦（耗竭型）n=%d = %.1f%%；真正健康（高意识型）n=%d = %.1f%%\n",
              nrow(norm), n_pd, 100*n_pd/nrow(norm), n_hc, 100*n_hc/nrow(norm)))

  # ---- 7.3 潜在痛苦 vs 真正健康：多结局对比（🔴 2026-10-06 起落表 + 补 d 的 CI）----
  # 说明：本表的数值此前只 print 到日志、未落任何 CSV，故 E7 长期缺少可追溯表源；
  #       现落 T7_potential_distress_vs_healthy.csv，并把 d 的 bootstrap 95% CI 一并落盘。
  a <- norm %>% filter(potential_depression); b <- norm %>% filter(healthy_control)
  cat(sprintf("  两组对比（n=%d vs %d）:\n", nrow(a), nrow(b)))
  cmp <- function(v, lab) {
    x <- a[[v]]; y <- b[[v]]
    if (length(x) < 3 || length(y) < 3) return(NULL)
    tt <- t.test(x, y)
    s <- sqrt((var(x, na.rm=TRUE) + var(y, na.rm=TRUE)) / 2)
    dl <- cohen_d_ci(x, y)
    data.frame(cohort = nm, 指标 = lab,
               n_pd = length(x[!is.na(x)]), n_hc = length(y[!is.na(y)]),
               潜在痛苦 = round(mean(x, na.rm=TRUE), 2),
               真正健康 = round(mean(y, na.rm=TRUE), 2),
               Δ = round(mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE), 2),
               d = round((mean(x, na.rm=TRUE) - mean(y, na.rm=TRUE)) / s, 2),
               d_lo = round(unname(dl[1]), 2), d_hi = round(unname(dl[2]), 2),
               p = signif(tt$p.value, 3))
  }
  tt <- do.call(rbind, Filter(Negate(is.null), list(
    cmp("gad7_total", "GAD-7"), cmp("phq9_total", "PHQ-9"),
    cmp("H_score", "H 值［分组定义变量·循环，仅供核对］"), cmp("ibs_positive", "IBS 阳性"),
    cmp("chronic_pain", "慢性疼痛"), cmp("bmi", "BMI"), cmp("sleep_hour", "睡眠时长"))))
  print(tt, row.names = FALSE)
  PD[[nm]] <- tt

  # ---- 7.4 三组群（group_type）----
  gt <- d %>% filter(!is.na(group_type)) %>% group_by(group_type) %>%
    summarise(n = n(), phq9 = mean(phq9_total, na.rm=TRUE), gad7 = mean(gad7_total, na.rm=TRUE),
              H = mean(H_score, na.rm=TRUE),
              ibs = mean(ibs_positive, na.rm=TRUE)*100, .groups = "drop") %>%
    arrange(factor(group_type, levels = c("真正健康","潜在抑郁","典型抑郁")))
  cat("  三组群（group_type）:\n")
  print(as.data.frame(gt %>% mutate(across(c(phq9,gad7,H), ~round(.x,3)), ibs = round(ibs,1))))

  # ---- 7.5 增量效度（🔴 2026-10-06 起落表）：α₃ 在**单一症状量表 PHQ-9** 之上预测另一结局 ----
  # ⚠️ 与模块 22 的 T22_incremental_validity（对照 = 四金标准组合、仅 Prolific n=142）口径不同，
  #    两者是不同对照、不同样本，引用时必须同句并列，禁止合并成一句或互相比大小。
  inc_row <- function(dat, y, kind) {
    dd <- dat %>% drop_na(all_of(y), phq9_total, alpha3)
    n <- nrow(dd)
    if (n < 30) return(NULL)
    if (kind == "cont") {
      m0 <- lm(reformulate("phq9_total", y), data = dd)
      m1 <- lm(reformulate(c("phq9_total", "alpha3"), y), data = dd)
      dR2 <- summary(m1)$r.squared - summary(m0)$r.squared
      pv  <- anova(m0, m1)$`Pr(>F)`[2]
      ci  <- bs_delta(dd, y, "cont")
      b   <- summary(m1)$coefficients["alpha3", ]
      cat(sprintf("  %s：PHQ-9 单独 R²=%.4f → 加入 α₃ 后 R²=%.4f（ΔR²=%+.4f，95%%CI %+.4f–%+.4f，F p=%s，n=%d）\n",
                  y, summary(m0)$r.squared, summary(m1)$r.squared, dR2,
                  unname(ci[1]), unname(ci[2]), signif(pv, 3), n))
      data.frame(cohort = nm, outcome = y, type = "连续（分层回归）", n = n,
                 R2_base = round(summary(m0)$r.squared, 4), R2_full = round(summary(m1)$r.squared, 4),
                 delta_R2 = round(dR2, 4),
                 delta_R2_lo = round(unname(ci[1]), 4), delta_R2_hi = round(unname(ci[2]), 4),
                 delta_R2_p = signif(pv, 3),
                 alpha_coef = round(b[1], 4), alpha_p = signif(b[4], 3),
                 AUC_base = NA_real_, AUC_full = NA_real_, delta_AUC = NA_real_,
                 delta_AUC_lo = NA_real_, delta_AUC_hi = NA_real_, delta_AUC_p = NA_real_)
    } else {
      m0 <- suppressWarnings(glm(reformulate("phq9_total", y), data = dd, family = binomial))
      m1 <- suppressWarnings(glm(reformulate(c("phq9_total", "alpha3"), y), data = dd, family = binomial))
      mn <- suppressWarnings(glm(reformulate("1", y), data = dd, family = binomial))
      lln <- as.numeric(logLik(mn))
      r2_0 <- 1 - as.numeric(logLik(m0))/lln; r2_1 <- 1 - as.numeric(logLik(m1))/lln
      a0 <- as.numeric(pROC::auc(pROC::roc(dd[[y]], predict(m0, type = "response"), quiet = TRUE)))
      a1 <- as.numeric(pROC::auc(pROC::roc(dd[[y]], predict(m1, type = "response"), quiet = TRUE)))
      pv <- anova(m0, m1, test = "Chisq")$`Pr(>Chi)`[2]
      ci <- bs_delta(dd, y, "bin")
      b  <- summary(m1)$coefficients["alpha3", ]
      cat(sprintf("  %s：PHQ-9 单独 AUC=%.4f → 加入 α₃ 后 AUC=%.4f（ΔAUC=%+.4f，95%%CI %+.4f–%+.4f，Chisq p=%s，n=%d）\n",
                  y, a0, a1, a1 - a0, unname(ci[1]), unname(ci[2]), signif(pv, 3), n))
      data.frame(cohort = nm, outcome = y, type = "二分类（ΔAUC）", n = n,
                 R2_base = round(r2_0, 4), R2_full = round(r2_1, 4), delta_R2 = round(r2_1 - r2_0, 4),
                 delta_R2_lo = NA_real_, delta_R2_hi = NA_real_, delta_R2_p = NA_real_,
                 alpha_coef = round(b[1], 4), alpha_p = signif(b[4], 3),
                 AUC_base = round(a0, 4), AUC_full = round(a1, 4), delta_AUC = round(a1 - a0, 4),
                 delta_AUC_lo = round(unname(ci[1]), 4), delta_AUC_hi = round(unname(ci[2]), 4),
                 delta_AUC_p = signif(pv, 3))
    }
  }
  IV[[nm]] <- do.call(rbind, Filter(Negate(is.null),
                                    list(inc_row(d, "gad7_total", "cont"),
                                         inc_row(d, "ibs_positive", "bin"))))

  # ---- 7.6 双正常亚组（🔴 2026-10-06 新增）：PHQ-9<10 且 GAD-7<10 ----
  # 目的：把「阈下识别」推进到**两个症状量表均低于临床阈值**的更严格子集。
  # 🔴 实测结果与预设不符，必须如实报告（不得粉饰、不得事后改标签）：
  #    组间 PHQ-9 / GAD-7 差**并未压到 0**（d ≈ 1.0–1.9）→
  #    ① 「结构分层与症状解耦」在本口径下**不成立**；
  #    ② 故本表只支持「**阈下亚临床梯度**」的存在，
  #       **不得**写成「测到完全无症状者的隐藏问题」。
  # 另：H 值是**分组定义变量**，其组间差属循环量
  #    同源比较仅供核对，不得作为独立外部效度结果引用。
  dn <- d %>% filter(phq9_total < 10, gad7_total < 10)
  a2 <- dn %>% filter(potential_depression); b2 <- dn %>% filter(healthy_control)
  cat(sprintf("  双正常者（PHQ-9<10 且 GAD-7<10）n=%d；潜在痛苦 n=%d / 真正健康 n=%d\n",
              nrow(dn), nrow(a2), nrow(b2)))
  cmp2 <- function(v, lab) {
    x <- a2[[v]]; y <- b2[[v]]
    if (length(x[!is.na(x)]) < 3 || length(y[!is.na(y)]) < 3) return(NULL)
    tt <- t.test(x, y)
    s <- sqrt((var(x, na.rm = TRUE) + var(y, na.rm = TRUE)) / 2)
    dl <- cohen_d_ci(x, y)
    data.frame(cohort = nm, 指标 = lab,
               n_pd = length(x[!is.na(x)]), n_hc = length(y[!is.na(y)]),
               潜在痛苦 = round(mean(x, na.rm = TRUE), 2),
               真正健康 = round(mean(y, na.rm = TRUE), 2),
               Δ = round(mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE), 2),
               d = round((mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)) / s, 2),
               d_lo = round(unname(dl[1]), 2), d_hi = round(unname(dl[2]), 2),
               p = signif(tt$p.value, 3))
  }
  tt2 <- do.call(rbind, Filter(Negate(is.null), list(
    cmp2("phq9_total", "PHQ-9［实测·非解耦］"), cmp2("gad7_total", "GAD-7［实测·非解耦］"),
    cmp2("H_score", "H 值［分组定义变量·循环，仅供核对］"), cmp2("ibs_positive", "IBS 阳性"),
    cmp2("chronic_pain", "慢性疼痛"), cmp2("bmi", "BMI"), cmp2("sleep_hour", "睡眠时长"))))
  print(tt2, row.names = FALSE)
  # 🔴 解耦与否的自证判读（如实打印，避免正文误读本表）
  if (!is.null(tt2)) {
    rr <- function(lab) tt2$d[tt2$指标 == lab]
    dp <- rr("PHQ-9［实测·非解耦］"); dg <- rr("GAD-7［实测·非解耦］")
    if (length(dp) && length(dg)) {
      cat(sprintf("    解耦自证：|d_PHQ9| = %.2f ｜ |d_GAD7| = %.2f → %s\n",
                  abs(dp), abs(dg),
                  if (abs(dp) < 0.2 && abs(dg) < 0.2)
                    "两组症状近似持平，支持解耦" else
                    "🔴 两组症状仍显著不同 → 「结构与症状解耦」不成立，本表只支持阈下亚临床梯度"))
    }
  }
  DN[[nm]] <- tt2

  OUT[[nm]] <- g %>% mutate(cohort = nm, .before = 1)
  GT[[nm]]  <- gt %>% mutate(cohort = nm, .before = 1)
}
save_tab(do.call(rbind, OUT), "T7_alpha_group")
save_tab(do.call(rbind, GT),  "T7_group_type")
save_tab(do.call(rbind, PD),  "T7_potential_distress_vs_healthy")
save_tab(do.call(rbind, IV),  "T7_incremental_vs_phq9")
save_tab(do.call(rbind, DN),  "T7_doubly_normal_split")
