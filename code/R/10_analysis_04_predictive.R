# =============================================================
# 10_analysis_04_predictive.R —— R3 预测效度：症状与疾病的 AUC / OR / R²
# 仅评价后台输出读数与外部结局的统计关联，不检查内部生成过程
# 两队列均适用
# =============================================================
# --- Windows 中文路径守卫：R 默认 locale 为 C 时无法解析中文路径 ---
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({a <- grep("^--file=", commandArgs(FALSE), value=TRUE); if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork=FALSE)) else getwd()})
source(file.path(P3_SCRIPT_DIR, "10_analysis_00_common.R"))
L <- load_cohorts(); AC <- c("alpha1","alpha2","alpha3","alpha4")
rows <- list()
for (nm in names(L)) {
  d <- L[[nm]]
  # 连续结局：α 联合预测的解释力
  for (out in c("phq9_total","gad7_total")) {
    dd <- d %>% select(all_of(c(out, AC))) %>% drop_na()
    f <- lm(as.formula(paste(out, "~", paste(AC, collapse="+"))), data = dd)
    r2 <- summary(f)$r.squared
    cat(sprintf("  %s: %s ~ 4α  R²=%.3f  (n=%d)\n", nm, out, r2, nrow(dd)))
    rows[[length(rows)+1]] <- data.frame(cohort=nm, outcome=out, type="R2", est=r2, n=nrow(dd))
  }
  # 二分结局：单 α AUC
  for (oc in c("ibs_positive")) {
    for (a in AC) {
      v <- auc_oriented(d[[oc]], d[[a]])
      rows[[length(rows)+1]] <- data.frame(cohort=nm, outcome=oc, type=paste0("AUC_",a),
                                           est=unname(v["auc"]), n=unname(v["n"]))
      cat(sprintf("  %s: %s ← %s  AUC=%.3f\n", nm, oc, a, v["auc"]))
    }
    cv <- cv_auc(d, AC, ycol = oc)
    rows[[length(rows)+1]] <- data.frame(cohort=nm, outcome=oc, type="AUC_4alpha_CV",
                                         est=cv["mean"], n=cv["n"])
    cat(sprintf("  %s: %s ← 4α 联合 5折CV  AUC=%.3f±%.3f  (n=%.0f)\n", nm, oc, cv["mean"], cv["sd"], cv["n"]))
  }
  # 疾病 OR（每 +0.1 α₃）
  for (dz in c("ibs_positive")) {
    dd <- d %>% mutate(a3_10 = alpha3/0.1) %>% select(all_of(c(dz,"a3_10"))) %>% drop_na()
    m <- suppressWarnings(glm(as.formula(paste(dz, "~ a3_10")), data = dd, family = binomial))
    b <- summary(m)$coefficients
    cat(sprintf("  %s: %s ← α₃ 每 +0.10  OR=%.2f [%.2f, %.2f]  p=%.3f\n", nm, dz,
                exp(b[2,1]), exp(b[2,1]-1.96*b[2,2]), exp(b[2,1]+1.96*b[2,2]), b[2,4]))
    rows[[length(rows)+1]] <- data.frame(cohort=nm, outcome=paste0(dz,"_OR_per0.1a3"),
                                         type="OR", est=exp(b[2,1]), n=nrow(dd))
  }
}
save_tab(do.call(rbind, rows), "T3_predictive_validity")
