# =============================================================
# 10_analysis_34_disease_or.R —— 七疾病 Logistic（α1–α4 预测自报慢病）
# 对应研究方案 A2（预测效度·七疾病 OR）/ C11（α₃ 与慢性疼痛负相关）/
#   C15（α₂ 与四类自报慢性病负相关）。此前管线仅 IBS(module04) 与
#   肥胖(module25) 有 OR 表，其余五病（慢性疼痛/高血压/糖尿病/高血脂/心脏病）
#   缺统一产出表 —— 本模块补齐。
# 两队列均适用；落盘 T34_disease_or.csv
# =============================================================
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))

ALPHAS  <- c("alpha1","alpha2","alpha3","alpha4")
COVARS  <- c("age_group","gender","education")
DIS7    <- c("ibs_positive","chronic_pain","heart_disease","hypertension",
            "diabetes","hyperlipidemia","obese")
DIS7_LAB<- c("IBS","慢性疼痛","心脏病","高血压","糖尿病","高血脂","肥胖")

L <- if (P3_COHORT %in% c("prolific","jianshu")) load_cohorts(P3_COHORT) else load_both()
rows <- list()

fit_or <- function(d, y, x, covars=NULL, unit, lab) {
  dd <- d %>% select(all_of(c(y, x, covars))) %>% drop_na()
  if (nrow(dd) < 40 || length(unique(dd[[y]])) < 2) return(NULL)
  f <- as.formula(paste(y, "~", paste(c(x, covars), collapse=" + ")))
  m <- suppressWarnings(glm(f, data=dd, family=binomial))
  b <- coef(m)[x]
  ci <- tryCatch(suppressMessages(confint(m))[x, ], error=function(e) c(NA_real_,NA_real_))
  data.frame(cohort=nm, disease=y, disease_lab=DIS7_LAB[which(DIS7==y)],
            alpha=lab, model=if(is.null(covars)) "crude" else "adjusted",
            unit=unit, OR=round(exp(b*unit),4),
            lo=round(exp(ci[1]*unit),4), hi=round(exp(ci[2]*unit),4),
            p=format.pval(summary(m)$coefficients[x,4], digits=3),
            n=nrow(dd), event_rate=round(mean(dd[[y]], na.rm=TRUE),4),
            stringsAsFactors=FALSE)
}

for (nm in names(L)) {
  d <- L[[nm]]
  cat(sprintf("\n=== %s：七疾病 OR（α1–α4）===\n", nm))
  for (dz in DIS7) {
    for (a in ALPHAS) {
      sdv <- sd(d[[a]], na.rm=TRUE)
      for (u in list(c("per+0.10",0.10), c("per+1SD",sdv))) {
        r1 <- fit_or(d, dz, a, NULL, as.numeric(u[2]), a)              # crude
        r2 <- fit_or(d, dz, a, COVARS, as.numeric(u[2]), a)           # adjusted
        if (!is.null(r1)) rows[[length(rows)+1]] <- r1
        if (!is.null(r2)) rows[[length(rows)+1]] <- r2
      }
    }
    cat(sprintf("  %s：事件率=%.3f 已完成 4α×2口径×2单位\n", DIS7_LAB[which(DIS7==dz)], mean(d[[dz]],na.rm=TRUE)))
  }
}
out <- do.call(rbind, rows)
save_tab(out, "T34_disease_or")
cat(sprintf("\n=== 34 模块完成：T34_disease_or.csv（%d 行）已写出 ===\n", nrow(out)))
