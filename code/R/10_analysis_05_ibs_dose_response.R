# =============================================================
# 10_analysis_05_ibs_dose_response.R —— R4 α₃-IBS 剂量反应
#   表1：低/高 α₃ 组间卡方 + 风险比(RR)
#   表2：α₃ 四分位 Logistic（OR，以 Q1 为参照；重点 Q4 vs Q1）
#   表3：α₃ 三分位序位趋势 + 连续 α₃ OR + Q1 vs Q3 率差
# 对应研究方案 A8 / C10（α₃ 与 IBS 负相关）+ 用户指定的卡方/趋势/OR 全套
# 两队列均适用；落盘 T5_ibs_dose_response.csv
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts()
all_rows <- list()

for (nm in names(L)) {
  d <- L[[nm]] %>% filter(is.finite(alpha3), is.finite(ibs_positive))
  d <- d %>% mutate(t3 = ntile(alpha3, 3), q4 = ntile(alpha3, 4))
  cat(sprintf("\n=== %s：α₃-IBS 剂量反应（n=%d）===\n", nm, nrow(d)))

  g <- d %>% group_by(t3) %>% summarise(n = n(), ibs = sum(ibs_positive, na.rm=TRUE),
                                        rate = mean(ibs_positive, na.rm=TRUE)*100,
                                        a3 = mean(alpha3, na.rm=TRUE),
                                        ibs_n = sum(ibs_positive, na.rm=TRUE),
                                        nonibs_n = n() - sum(ibs_positive, na.rm=TRUE),
                                        .groups="drop")
  print(as.data.frame(g))

  # ---- 表3a：趋势检验（三分位序数）----
  m <- suppressWarnings(glm(ibs_positive ~ t3, data = d, family = binomial))
  b <- summary(m)$coefficients
  or_trend <- exp(b[2,1]); p_trend <- b[2,4]
  cat(sprintf("  趋势检验（三分位序数）：OR=%.2f  p=%.4f\n", or_trend, p_trend))

  # ---- 表3b：连续 α₃ ----
  m2 <- suppressWarnings(glm(ibs_positive ~ alpha3, data = d, family = binomial))
  b2 <- summary(m2)$coefficients
  or_cont <- exp(b2[2,1]/10); p_cont <- b2[2,4]
  cat(sprintf("  连续 α₃：OR(每+0.1)=%.2f  p=%.4f\n", or_cont, p_cont))

  # ---- 表3c：Q1 vs Q3 阳性率差 ----
  q13 <- g$rate[1] - g$rate[3]
  cat(sprintf("  Q1 vs Q3 阳性率差 = %.1f 个百分点\n", q13))

  # ---- 表1：低/高 α₃ 二分（中位切分）卡方 + 风险比 ----
  med <- median(d$alpha3, na.rm=TRUE)
  d$bin <- ifelse(d$alpha3 < med, "低α₃", "高α₃")
  tab <- table(d$bin, d$ibs_positive)
  chi <- suppressWarnings(chisq.test(tab))
  rate_low  <- mean(d$ibs_positive[d$bin=="低α₃"], na.rm=TRUE)
  rate_high <- mean(d$ibs_positive[d$bin=="高α₃"], na.rm=TRUE)
  rr <- if (rate_low>0) rate_high/rate_low else NA_real_
  cat(sprintf("  低/高 α₃ 卡方：χ²=%.2f  p=%.4f ｜ 低α₃阳性率=%.1f%% 高α₃阳性率=%.1f%% ｜ RR(高/低)=%.2f\n",
              chi$statistic, chi$p.value, rate_low*100, rate_high*100, rr))

  # ---- 表2：四分位 Logistic（Q1 参照），重点 Q4 vs Q1 ----
  d$q4f <- factor(d$q4, levels=1:4, labels=c("Q1","Q2","Q3","Q4"))
  mq <- suppressWarnings(glm(ibs_positive ~ q4f, data = d, family = binomial))
  cmq <- summary(mq)$coefficients
  ciq <- tryCatch(suppressMessages(confint(mq)), error=function(e) matrix(NA, nrow=nrow(cmq), ncol=2))
  # 各分位 vs Q1 的 OR（统一 15 列 schema，避免 rbind 列数不一致）
  ors <- exp(c(cmq["q4fQ2",1],cmq["q4fQ3",1],cmq["q4fQ4",1]))
  los <- exp(c(ciq["q4fQ2",1],ciq["q4fQ3",1],ciq["q4fQ4",1]))
  his <- exp(c(ciq["q4fQ2",2],ciq["q4fQ3",2],ciq["q4fQ4",2]))
  ps  <- format.pval(c(cmq["q4fQ2",4],cmq["q4fQ3",4],cmq["q4fQ4",4]), digits=3)
  qOR <- data.frame(
    cohort=nm, analysis="表2_四分位OR_vs_Q1", cut="四分位(Q1参照)",
    group1="Q1", n1=sum(d$q4==1), rate1_pct=round(mean(d$ibs_positive[d$q4==1])*100,2),
    group2=c("Q2","Q3","Q4"), n2=sapply(2:4,function(i) sum(d$q4==i)),
    rate2_pct=sapply(2:4,function(i) round(mean(d$ibs_positive[d$q4==i])*100,2)),
    stat=NA_real_, p=ps,
    OR=round(ors,3), lo=round(los,3), hi=round(his,3),
    extra=c("OR(Q2/Q1)","OR(Q3/Q1)","OR(Q4/Q1)"), stringsAsFactors=FALSE)
  # 单独报 OR(Q4 vs Q1)（用户指定项）
  or_q4 <- exp(cmq["q4fQ4",1]); lo_q4 <- exp(ciq["q4fQ4",1]); hi_q4 <- exp(ciq["q4fQ4",2]); p_q4 <- cmq["q4fQ4",4]
  cat(sprintf("  OR(Q4 vs Q1)=%.2f [%.2f, %.2f]  p=%.4f\n", or_q4, lo_q4, hi_q4, p_q4))

  # ---- 组装落盘表 ----
  all_rows[[length(all_rows)+1]] <- data.frame(
    cohort=nm, analysis="表1_低高α3卡方_RR", cut=sprintf("中位=%.3f", med),
    group1="低α₃", n1=sum(d$bin=="低α₃"), rate1_pct=round(rate_low*100,2),
    group2="高α₃", n2=sum(d$bin=="高α₃"), rate2_pct=round(rate_high*100,2),
    stat=round(chi$statistic,3), p=round(chi$p.value,4),
    OR=round(rr,3), lo=NA_real_, hi=NA_real_, extra="RR(高/低)", stringsAsFactors=FALSE)
  all_rows[[length(all_rows)+1]] <- data.frame(
    cohort=nm, analysis="表2_OR_Q4vsQ1", cut="四分位(Q1参照)",
    group1="Q1", n1=sum(d$q4==1), rate1_pct=round(mean(d$ibs_positive[d$q4==1],na.rm=TRUE)*100,2),
    group2="Q4", n2=sum(d$q4==4), rate2_pct=round(mean(d$ibs_positive[d$q4==4],na.rm=TRUE)*100,2),
    stat=NA_real_, p=round(p_q4,4),
    OR=round(or_q4,3), lo=round(lo_q4,3), hi=round(hi_q4,3), extra="OR(Q4/Q1)", stringsAsFactors=FALSE)
  all_rows[[length(all_rows)+1]] <- qOR
  all_rows[[length(all_rows)+1]] <- data.frame(
    cohort=nm, analysis="表3_趋势_三分位序数", cut="",
    group1="", n1=NA_integer_, rate1_pct=NA_real_,
    group2="", n2=NA_integer_, rate2_pct=NA_real_,
    stat=NA_real_, p=round(p_trend,4), OR=round(or_trend,3), lo=NA_real_, hi=NA_real_,
    extra="趋势OR/每增1分位", stringsAsFactors=FALSE)
  all_rows[[length(all_rows)+1]] <- data.frame(
    cohort=nm, analysis="表3_连续α3_OR", cut="per+0.1",
    group1="", n1=NA_integer_, rate1_pct=NA_real_,
    group2="", n2=NA_integer_, rate2_pct=NA_real_,
    stat=NA_real_, p=round(p_cont,4), OR=round(or_cont,3), lo=NA_real_, hi=NA_real_,
    extra="OR/每+0.1α3", stringsAsFactors=FALSE)
  all_rows[[length(all_rows)+1]] <- data.frame(
    cohort=nm, analysis="表3_Q1vsQ3率差", cut="",
    group1="Q1", n1=sum(d$t3==1), rate1_pct=round(g$rate[1],2),
    group2="Q3", n2=sum(d$t3==3), rate2_pct=round(g$rate[3],2),
    stat=NA_real_, p=NA_real_, OR=NA_real_, lo=NA_real_, hi=NA_real_,
    extra=sprintf("率差=%.1fpp", q13), stringsAsFactors=FALSE)
}
save_tab(do.call(rbind, all_rows), "T5_ibs_dose_response")
cat("\n=== 05 模块完成：T5_ibs_dose_response.csv 已写出 ===\n")
