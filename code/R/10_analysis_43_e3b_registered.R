# E3b registered C1/C2/C3, output-level public analysis only.
# C1 uses the registered Cliff's-delta margin, not a raw mean margin.
suppressWarnings(Sys.setlocale("LC_ALL", "zh_CN.UTF-8"))
P3_SCRIPT_DIR <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]), mustWork = FALSE)) else getwd()
})
PROJ_E3B <- normalizePath(Sys.getenv("P3_ROOT", unset = file.path(P3_SCRIPT_DIR, "..", "..")),
                          winslash = "/", mustWork = FALSE)
DERIVED_E3B <- file.path(PROJ_E3B, "code", "data", "derived")
OUT_E3B <- Sys.getenv("P3_E3B_OUT", unset = "")
if (!nzchar(OUT_E3B)) {
  parent_out <- Sys.getenv("P3_OUT_DIR", unset = "")
  OUT_E3B <- if (nzchar(parent_out)) file.path(parent_out, "E3b_results") else
    file.path(PROJ_E3B, "data_output", "E3b_results")
}
dir.create(OUT_E3B, recursive = TRUE, showWarnings = FALSE)

read_public_e3b <- function(name) {
  read.csv(file.path(DERIVED_E3B, name), stringsAsFactors = FALSE,
           check.names = FALSE, fileEncoding = "UTF-8-BOM", na.strings = c("", "NA"))
}
save_e3b <- function(x, name) {
  write.csv(x, file.path(OUT_E3B, name), row.names = FALSE,
            fileEncoding = "UTF-8", na = "")
  cat("saved", file.path(OUT_E3B, name), "\n")
}

alpha <- paste0("alpha", 1:4)
long <- read_public_e3b("p3_base_phase2_long.csv")
reference <- read_public_e3b("p3_base_jianshu.csv")
for (c in c(alpha, "phq9_total", "ibs_positive")) {
  if (c %in% names(long)) long[[c]] <- as.numeric(long[[c]])
  if (c %in% names(reference)) reference[[c]] <- as.numeric(reference[[c]])
}
long$batch <- as.integer(long$batch)
b2 <- long[long$batch == 2 & long$questionnaire == "新20", , drop = FALSE]
stopifnot(nrow(b2) == 240, nrow(reference) == 504)

# Direct port of statsmodels' rank_compare_2indep midrank/Brunner-Munzel
# moments, including its Welch-Satterthwaite degrees of freedom.
rank_compare_e3b <- function(x1, x2) {
  n1 <- length(x1); n2 <- length(x2); n <- n1 + n2
  r <- rank(c(x1, x2), ties.method = "average")
  r1 <- r[seq_len(n1)]; r2 <- r[n1 + seq_len(n2)]
  ri1 <- rank(x1, ties.method = "average")
  ri2 <- rank(x2, ties.method = "average")
  mr1 <- mean(r1); mr2 <- mean(r2)
  mir1 <- mean(ri1); mir2 <- mean(ri2)
  s1 <- sum((r1 - ri1 - mr1 + mir1)^2) / (n1 - 1)
  s2 <- sum((r2 - ri2 - mr2 + mir2)^2) / (n2 - 1)
  numerator <- (n1 * s1 + n2 * s2)^2
  denominator <- (n1 * s1)^2 / (n1 - 1) + (n2 * s2)^2 / (n2 - 1)
  df <- numerator / denominator
  var1 <- s1 / (n - n1)^2
  var2 <- s2 / (n - n2)^2
  var_prob <- var1 / n1 + var2 / n2
  prob1 <- (mr1 - (n1 + 1) / 2) / n2
  list(prob1 = prob1, var_prob = var_prob, df = df)
}
rank_tost_e3b <- function(a, b, margin = 0.147) {
  fit <- rank_compare_e3b(a, b)
  delta <- 2 * fit$prob1 - 1
  se <- 2 * sqrt(fit$var_prob)
  crit <- qt(0.95, fit$df)
  lo <- delta - crit * se; hi <- delta + crit * se
  p_lower <- pt((delta + margin) / se, fit$df, lower.tail = FALSE)
  p_upper <- pt((delta - margin) / se, fit$df, lower.tail = TRUE)
  c(delta = delta, se = se, lo = lo, hi = hi, df = fit$df,
    p_lower = p_lower, p_upper = p_upper,
    equivalent = max(p_lower, p_upper) < 0.05)
}

c1_rows <- list(); distribution_rows <- list(); idx <- 1L; didx <- 1L
for (variable in alpha) {
  a <- b2[[variable]]; a <- a[is.finite(a)]
  b <- reference[[variable]]; b <- b[is.finite(b)]
  z <- rank_tost_e3b(a, b)
  c1_rows[[idx]] <- data.frame(
    metric = "C1_registered_rank_TOST",
    group_a = "second_wave_batch2_new20",
    group_b = "existing_Chinese_cohort_first20",
    variable = variable, n_a = length(a), n_b = length(b),
    cliff_delta = unname(z["delta"]), se = unname(z["se"]), df = unname(z["df"]),
    ci90_lo = unname(z["lo"]), ci90_hi = unname(z["hi"]),
    p_tost_lower = unname(z["p_lower"]), p_tost_upper = unname(z["p_upper"]),
    margin = 0.147, equivalent = as.logical(z["equivalent"]),
    analysis_role = "registered_confirmatory_C1",
    method = "Brunner-Munzel rank inference; delta=2P(a>b with midrank ties)-1",
    stringsAsFactors = FALSE)
  idx <- idx + 1L
  for (pair in list(list("second_wave_batch2_new20", a),
                    list("existing_Chinese_cohort_first20", b))) {
    label <- pair[[1]]; values <- pair[[2]]
    distribution_rows[[didx]] <- data.frame(
      variable = variable, group = label, n = length(values),
      mean = mean(values), sd = sd(values), min = min(values),
      p25 = unname(quantile(values, 0.25)), median = median(values),
      p75 = unname(quantile(values, 0.75)), max = max(values),
      analysis_role = "registered_distribution_description", stringsAsFactors = FALSE)
    didx <- didx + 1L
  }
}
c1 <- do.call(rbind, c1_rows)
save_e3b(c1, "C1_预注册标准化主结果.csv")
save_e3b(do.call(rbind, distribution_rows), "C1_分布描述.csv")

c2_rows <- list(); idx <- 1L
samples <- list(
  list("existing_Chinese_cohort", reference, "direction_reference"),
  list("second_wave_batch2", b2, "registered_directional_replication"))
for (entry in samples) {
  label <- entry[[1]]; frame <- entry[[2]]; role <- entry[[3]]
  phq <- frame[is.finite(frame$alpha3) & is.finite(frame$phq9_total), ]
  phq_test <- suppressWarnings(cor.test(phq$alpha3, phq$phq9_total,
                                        method = "spearman", exact = FALSE))
  c2_rows[[idx]] <- data.frame(
    sample = label, outcome = "PHQ-9_total", n = nrow(phq),
    statistic = "Spearman_rho_alpha3_outcome",
    estimate = unname(phq_test$estimate), p_two_sided = phq_test$p.value,
    expected_direction = "negative", direction_replicated = phq_test$estimate < 0,
    analysis_role = role, n_positive = NA_integer_, stringsAsFactors = FALSE)
  idx <- idx + 1L
  ibs <- frame[is.finite(frame$alpha3) & is.finite(frame$ibs_positive), ]
  y <- as.integer(ibs$ibs_positive)
  ibs_test <- cor.test(ibs$alpha3, y, method = "pearson")
  c2_rows[[idx]] <- data.frame(
    sample = label, outcome = "IBS_positive", n = nrow(ibs),
    statistic = "point_biserial_r_alpha3_outcome",
    estimate = unname(ibs_test$estimate), p_two_sided = ibs_test$p.value,
    expected_direction = "negative", direction_replicated = ibs_test$estimate < 0,
    analysis_role = role, n_positive = sum(y), stringsAsFactors = FALSE)
  idx <- idx + 1L
}
c2 <- do.call(rbind, c2_rows)
save_e3b(c2, "C2_方向复现.csv")
c2_replication <- c2[c2$analysis_role == "registered_directional_replication", ]
stopifnot(nrow(c2_replication) == 2)

c3 <- read_public_e3b("p3_c3_option_distribution.csv")
stopifnot(nrow(c3) == 100, length(unique(c3$item_id)) == 20,
          identical(sort(unique(as.integer(c3$option_position))), 1:5),
          all(tapply(c3$n, c3$item_id, sum) == 240),
          all(abs(c3$percent - 100 * c3$n / c3$source_n) < 1e-5))
save_e3b(c3, "C3_选项级分布.csv")

summary_rows <- list(); idx <- 1L
for (i in seq_len(nrow(c1))) {
  row <- c1[i, ]
  summary_rows[[idx]] <- data.frame(
    analysis = "C1", registration_status = "registered",
    group_a = row$group_a, group_b = row$group_b, variable = row$variable,
    estimand = "Cliff_delta", estimate = row$cliff_delta,
    ci90_lo = row$ci90_lo, ci90_hi = row$ci90_hi,
    p_tost_lower = row$p_tost_lower, p_tost_upper = row$p_tost_upper,
    criterion_met = row$equivalent, stringsAsFactors = FALSE)
  idx <- idx + 1L
}
for (i in seq_len(nrow(c2_replication))) {
  row <- c2_replication[i, ]
  summary_rows[[idx]] <- data.frame(
    analysis = "C2", registration_status = "registered_direction",
    group_a = row$sample, group_b = row$outcome, variable = "alpha3",
    estimand = row$statistic, estimate = row$estimate,
    ci90_lo = NA_real_, ci90_hi = NA_real_, p_tost_lower = NA_real_,
    p_tost_upper = NA_real_, criterion_met = row$direction_replicated,
    stringsAsFactors = FALSE)
  idx <- idx + 1L
}
save_e3b(do.call(rbind, summary_rows), "E3b_主结果.csv")

c1_passed <- c1$variable[c1$equivalent]
c1_failed <- c1$variable[!c1$equivalent]
c2_pass <- all(c2_replication$direction_replicated)
summary <- c("# E3b registered-analysis summary", "",
  "C1 compares second-wave batch 2 (n=240) with the existing Chinese cohort (n=504).",
  "The margin is |Cliff's delta| < 0.147; it is not a raw mean-difference bound.", "",
  "| alpha | Cliff delta | 90% CI | lower p | upper p | equivalent |",
  "|---|---:|---|---:|---:|---|")
for (i in seq_len(nrow(c1))) {
  summary <- c(summary, sprintf("| %s | %.6f | [%.6f, %.6f] | %.6g | %.6g | %s |",
    c1$variable[i], c1$cliff_delta[i], c1$ci90_lo[i], c1$ci90_hi[i],
    c1$p_tost_lower[i], c1$p_tost_upper[i], if (isTRUE(c1$equivalent[i])) "TRUE" else "FALSE"))
}
summary <- c(summary, "",
  sprintf("C1 is partially supported: passed=%s; not passed=%s.",
          paste(c1_passed, collapse = ","), paste(c1_failed, collapse = ",")),
  sprintf("C2 direction-only descriptive replication (both expected negative): direction_consistent=%s; no registered model or significance threshold.", if (isTRUE(c2_pass)) "TRUE" else "FALSE"),
  "C3 comprises 20 items x 5 option positions, 4,800 responses; minimum cell n=2; descriptive only, with no registered sparsity failure threshold.",
  "The batch-1 versus batch-2 mean-scale table is post hoc and is not registered C1.",
  "E3a and order analyses are protocol-specified but unregistered and non-confirmatory.",
  "No item text, mapping, scoring key, parameter, or backend formula is disclosed.")
writeLines(summary, file.path(OUT_E3B, "E3b_分析摘要.md"), useBytes = TRUE)
cat("C1 passed:", paste(c1_passed, collapse = ","), "; failed:",
    paste(c1_failed, collapse = ","), "; C2 directions passed:", c2_pass, "\n")
