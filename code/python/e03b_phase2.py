# -*- coding: utf-8 -*-
"""E3b/E3a downstream analyses from public output-level data.

Registered E3b (2026-10-10):
  C1: second-wave batch 2 (new item set, n=240) versus the existing
      Chinese-language cohort (first item set, n=504). The confirmatory
      equivalence margin is |Cliff's delta| < 0.147, evaluated by rank-based
      Brunner-Munzel inference and two one-sided tests. The value 0.147 is
      not a raw mean-difference margin.
  C2: directional replication of alpha3--IBS and alpha3--PHQ-9 associations.
  C3: descriptive option-position distributions for the 20 new items.

The earlier batch-1 (n=166) versus batch-2 (n=240) mean-scale comparison is
post hoc and is never substituted for registered C1. E2/E3a and the order
check were specified in the approved protocol but not separately registered;
they remain non-confirmatory.

No item text, item-to-score mapping, proprietary parameter, backend formula,
or internal representation is read or disclosed by this module.
"""
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.nonparametric import rank_compare_2indep

ROOT = os.path.abspath(os.environ.get(
    "P3_ROOT", os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))))
DERIVED = os.path.join(ROOT, "code", "data", "derived")
DEFAULT_OUT = (os.path.join(os.environ["P3_RESULT_DIR"], "E3b_results")
               if os.environ.get("P3_RESULT_DIR")
               else os.path.join(ROOT, "data_output", "E3b_results"))
OUT = os.path.abspath(os.environ.get("P3_E3B_OUT", DEFAULT_OUT))
os.makedirs(OUT, exist_ok=True)

ALPHA = ["alpha1", "alpha2", "alpha3", "alpha4"]
CARR = ["carrier_a", "carrier_b", "carrier_c", "carrier_d"]
READINGS = ALPHA + CARR + ["H_score"]
MARGIN = 0.147


def read_public(name):
    return pd.read_csv(os.path.join(DERIVED, name), dtype=str,
                       keep_default_na=False, encoding="utf-8-sig")


long = read_public("p3_base_phase2_long.csv")
wide = read_public("p3_base_phase2_wide.csv")
reference = read_public("p3_base_jianshu.csv")
for frame, columns in ((long, READINGS + ["phq9_total", "ibs_positive"]),
                       (reference, READINGS + ["phq9_total", "ibs_positive"])):
    for col in columns:
        if col in frame:
            frame[col] = pd.to_numeric(frame[col], errors="coerce")
long["batch"] = pd.to_numeric(long["batch"], errors="raise").astype(int)
for col in READINGS:
    for suffix in ("_new", "_old"):
        if col + suffix in wide:
            wide[col + suffix] = pd.to_numeric(wide[col + suffix], errors="coerce")

b1 = long[(long.batch == 1) & (long.questionnaire == "新20")].copy()
b2 = long[(long.batch == 2) & (long.questionnaire == "新20")].copy()
assert len(b1) == 166 and len(b2) == 240 and len(reference) == 504
assert len(wide) == 166


def rank_tost(a, b, margin=MARGIN):
    """Rank-based TOST for Cliff's delta; sign is mean(a-b) direction."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    result = rank_compare_2indep(a, b, use_t=True)
    delta = 2 * result.prob1 - 1
    se = 2 * np.sqrt(result.var_prob)
    ci = result.conf_int(alpha=0.10)
    lo, hi = 2 * ci[0] - 1, 2 * ci[1] - 1
    df = result.df
    p_lower = stats.t.sf((delta + margin) / se, df)
    p_upper = stats.t.cdf((delta - margin) / se, df)
    equivalent = bool(max(p_lower, p_upper) < 0.05)
    return delta, se, lo, hi, df, p_lower, p_upper, equivalent


def pooled_sd(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    return np.sqrt(((len(a) - 1) * a.var(ddof=1) +
                    (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))


def mean_tost_independent(a, b, margin):
    """Welch mean-scale TOST; used only for explicitly post-hoc analyses."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    diff = a.mean() - b.mean()
    v1, v2 = a.var(ddof=1), b.var(ddof=1)
    se = np.sqrt(v1 / len(a) + v2 / len(b))
    df = (v1 / len(a) + v2 / len(b)) ** 2 / (
        (v1 / len(a)) ** 2 / (len(a) - 1) +
        (v2 / len(b)) ** 2 / (len(b) - 1))
    crit = stats.t.ppf(0.95, df)
    lo, hi = diff - crit * se, diff + crit * se
    p_lower = stats.t.sf((diff + margin) / se, df)
    p_upper = stats.t.cdf((diff - margin) / se, df)
    return diff, lo, hi, p_lower, p_upper, bool(max(p_lower, p_upper) < 0.05)


def paired_tost(diff, margin):
    diff = np.asarray(diff, float)
    mean = diff.mean()
    se = diff.std(ddof=1) / np.sqrt(len(diff))
    df = len(diff) - 1
    crit = stats.t.ppf(0.95, df)
    lo, hi = mean - crit * se, mean + crit * se
    p_lower = stats.t.sf((mean + margin) / se, df)
    p_upper = stats.t.cdf((mean - margin) / se, df)
    return mean, lo, hi, p_lower, p_upper, bool(max(p_lower, p_upper) < 0.05)


def icc21(data):
    """Two-way random-effects absolute-agreement ICC(2,1) for two forms."""
    data = np.asarray(data, float)
    n, k = data.shape
    grand = data.mean()
    ss_rows = k * np.sum((data.mean(axis=1) - grand) ** 2)
    ss_cols = n * np.sum((data.mean(axis=0) - grand) ** 2)
    ss_error = ((data - grand) ** 2).sum() - ss_rows - ss_cols
    msr = ss_rows / (n - 1)
    msc = ss_cols / (k - 1)
    mse = ss_error / ((n - 1) * (k - 1))
    return (msr - mse) / (msr + (k - 1) * mse + k * (msc - mse) / n)


def save(frame, name):
    path = os.path.join(OUT, name)
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    print("saved", path)


# C1: registered rank-scale comparison against the existing Chinese cohort.
c1_rows = []
distribution_rows = []
for variable in ALPHA:
    a = b2[variable].dropna().to_numpy()
    b = reference[variable].dropna().to_numpy()
    delta, se, lo, hi, df, pl, pu, equivalent = rank_tost(a, b)
    c1_rows.append({
        "metric": "C1_registered_rank_TOST",
        "group_a": "second_wave_batch2_new20",
        "group_b": "existing_Chinese_cohort_first20",
        "variable": variable, "n_a": len(a), "n_b": len(b),
        "cliff_delta": delta, "se": se, "df": df,
        "ci90_lo": lo, "ci90_hi": hi,
        "p_tost_lower": pl, "p_tost_upper": pu,
        "margin": MARGIN, "equivalent": equivalent,
        "analysis_role": "registered_confirmatory_C1",
        "method": "Brunner-Munzel rank inference; delta=2P(a>b with midrank ties)-1",
    })
    for label, values in (("second_wave_batch2_new20", a),
                          ("existing_Chinese_cohort_first20", b)):
        distribution_rows.append({
            "variable": variable, "group": label, "n": len(values),
            "mean": values.mean(), "sd": values.std(ddof=1),
            "min": values.min(), "p25": np.quantile(values, 0.25),
            "median": np.median(values), "p75": np.quantile(values, 0.75),
            "max": values.max(), "analysis_role": "registered_distribution_description",
        })
c1 = pd.DataFrame(c1_rows)
save(c1, "C1_预注册标准化主结果.csv")
save(pd.DataFrame(distribution_rows), "C1_分布描述.csv")


# C2: directional replication; no equivalence or clinical utility is implied.
c2_rows = []
for sample_label, frame, role in (
        ("existing_Chinese_cohort", reference, "direction_reference"),
        ("second_wave_batch2", b2, "registered_directional_replication")):
    phq = frame[["alpha3", "phq9_total"]].dropna()
    rho, phq_p = stats.spearmanr(phq.alpha3, phq.phq9_total)
    c2_rows.append({
        "sample": sample_label, "outcome": "PHQ-9_total", "n": len(phq),
        "statistic": "Spearman_rho_alpha3_outcome", "estimate": rho,
        "p_two_sided": phq_p, "expected_direction": "negative",
        "direction_replicated": bool(rho < 0), "analysis_role": role,
    })
    ibs = frame[["alpha3", "ibs_positive"]].dropna()
    y = ibs.ibs_positive.astype(int)
    rpb, ibs_p = stats.pointbiserialr(y, ibs.alpha3)
    c2_rows.append({
        "sample": sample_label, "outcome": "IBS_positive", "n": len(ibs),
        "n_positive": int(y.sum()), "statistic": "point_biserial_r_alpha3_outcome",
        "estimate": rpb, "p_two_sided": ibs_p, "expected_direction": "negative",
        "direction_replicated": bool(rpb < 0), "analysis_role": role,
    })
c2 = pd.DataFrame(c2_rows)
save(c2, "C2_方向复现.csv")
c2_replication = c2[c2.analysis_role == "registered_directional_replication"]
assert len(c2_replication) == 2


# C3: public aggregate counts only. Item wording and response keys are absent.
c3 = read_public("p3_c3_option_distribution.csv")
for col in ("option_position", "n", "percent", "source_n"):
    c3[col] = pd.to_numeric(c3[col], errors="raise")
assert len(c3) == 100 and c3.item_id.nunique() == 20
assert set(c3.option_position.unique()) == {1, 2, 3, 4, 5}
assert c3.groupby("item_id").n.sum().eq(240).all()
assert np.allclose(c3.percent, 100 * c3.n / c3.source_n, atol=1e-5)
save(c3, "C3_选项级分布.csv")


# Compact primary-result table for navigation.
summary_rows = []
for row in c1_rows:
    summary_rows.append({
        "analysis": "C1", "registration_status": "registered",
        "group_a": row["group_a"], "group_b": row["group_b"],
        "variable": row["variable"], "estimand": "Cliff_delta",
        "estimate": row["cliff_delta"], "ci90_lo": row["ci90_lo"],
        "ci90_hi": row["ci90_hi"], "p_tost_lower": row["p_tost_lower"],
        "p_tost_upper": row["p_tost_upper"], "criterion_met": row["equivalent"],
    })
for _, row in c2_replication.iterrows():
    summary_rows.append({
        "analysis": "C2", "registration_status": "registered_direction",
        "group_a": row["sample"], "group_b": row["outcome"],
        "variable": "alpha3", "estimand": row["statistic"],
        "estimate": row["estimate"], "ci90_lo": np.nan, "ci90_hi": np.nan,
        "p_tost_lower": np.nan, "p_tost_upper": np.nan,
        "criterion_met": row["direction_replicated"],
    })
save(pd.DataFrame(summary_rows), "E3b_主结果.csv")


# Post-hoc mean-scale comparison: batch 1 versus batch 2, not registered C1.
posthoc_rows = []
for variable in READINGS[:-1]:
    a = b1[variable].dropna().to_numpy()
    b = b2[variable].dropna().to_numpy()
    sp = pooled_sd(a, b)
    raw = mean_tost_independent(a, b, MARGIN)
    standardized = mean_tost_independent(a, b, 0.262 * sp)
    posthoc_rows.append({
        "variable": variable, "n_batch1": len(a), "n_batch2": len(b),
        "mean_batch1": a.mean(), "sd_batch1": a.std(ddof=1),
        "mean_batch2": b.mean(), "sd_batch2": b.std(ddof=1),
        "mean_diff_batch1_minus_batch2": raw[0], "pooled_sd": sp,
        "cohen_d": raw[0] / sp,
        "raw_margin": MARGIN, "raw_ci90_lo": raw[1], "raw_ci90_hi": raw[2],
        "raw_p_tost_lower": raw[3], "raw_p_tost_upper": raw[4],
        "raw_equivalent": raw[5], "standardized_margin": 0.262 * sp,
        "standardized_ci90_lo": standardized[1],
        "standardized_ci90_hi": standardized[2],
        "standardized_p_tost_lower": standardized[3],
        "standardized_p_tost_upper": standardized[4],
        "standardized_equivalent": standardized[5],
        "analysis_role": "post_hoc_batch1_vs_batch2_not_registered_C1",
        "boundary": "0.147 here is a raw mean-difference sensitivity margin, not Cliff's delta",
    })
save(pd.DataFrame(posthoc_rows), "PostHoc_第一批对第二批_均值TOST.csv")


# Protocol-specified but unregistered dual-form analysis.
e3a_rows = []
for variable in ALPHA:
    new = wide[variable + "_new"].to_numpy()
    old = wide[variable + "_old"].to_numpy()
    differences = new - old
    sd_diff = differences.std(ddof=1)
    raw = paired_tost(differences, MARGIN)
    standardized = paired_tost(differences, 0.262 * sd_diff)
    e3a_rows.append({
        "variable": variable, "n": len(differences),
        "mean_new": new.mean(), "sd_new": new.std(ddof=1),
        "mean_old": old.mean(), "sd_old": old.std(ddof=1),
        "raw_diff": differences.mean(), "sd_diff": sd_diff,
        "cohen_dz": differences.mean() / sd_diff,
        "pearson_r": np.corrcoef(new, old)[0, 1],
        "icc21_absolute_agreement": icc21(np.column_stack([new, old])),
        "raw_margin": MARGIN, "raw_ci90_lo": raw[1], "raw_ci90_hi": raw[2],
        "raw_p_tost_lower": raw[3], "raw_p_tost_upper": raw[4],
        "raw_equivalent": raw[5], "standardized_margin": 0.262 * sd_diff,
        "standardized_ci90_lo": standardized[1],
        "standardized_ci90_hi": standardized[2],
        "standardized_p_tost_lower": standardized[3],
        "standardized_p_tost_upper": standardized[4],
        "standardized_equivalent": standardized[5],
        "analysis_role": "protocol_specified_unregistered_nonconfirmatory_E3a",
    })
save(pd.DataFrame(e3a_rows), "E3a_配对标准化敏感性.csv")


# Unregistered order check under the historical raw mean-scale sensitivity rule.
order_rows = []
order_a = long[(long.batch == 1) & (long.order_group == "先新后旧") &
               (long.questionnaire == "新20")]
order_b = long[(long.batch == 1) & (long.order_group == "先旧后新") &
               (long.questionnaire == "新20")]
for variable in ALPHA:
    result = mean_tost_independent(order_a[variable], order_b[variable], MARGIN)
    order_rows.append({
        "variable": variable, "n_new_first": len(order_a), "n_old_first": len(order_b),
        "mean_diff": result[0], "ci90_lo": result[1], "ci90_hi": result[2],
        "p_tost_lower": result[3], "p_tost_upper": result[4],
        "equivalent_raw_0_147": result[5],
        "analysis_role": "protocol_specified_unregistered_order_check",
    })
save(pd.DataFrame(order_rows), "顺序效应_非注册TOST.csv")
save(b2[READINGS[:-1]].agg(["mean", "std", "min", "max"]).T.reset_index().rename(
    columns={"index": "variable"}), "D_第2批描述性.csv")


passed = c1.loc[c1.equivalent, "variable"].tolist()
failed = c1.loc[~c1.equivalent, "variable"].tolist()
c2_pass = c2_replication.direction_replicated.all()
summary = [
    "# E3b registered-analysis summary",
    "",
    "C1 compares second-wave batch 2 (n=240) with the existing Chinese cohort (n=504).",
    "The margin is |Cliff's delta| < 0.147; it is not a raw mean-difference bound.",
    "",
    "| alpha | Cliff delta | 90% CI | lower p | upper p | equivalent |",
    "|---|---:|---|---:|---:|---|",
]
for row in c1_rows:
    summary.append("| {variable} | {cliff_delta:.6f} | [{ci90_lo:.6f}, {ci90_hi:.6f}] | "
                   "{p_tost_lower:.6g} | {p_tost_upper:.6g} | {equivalent} |".format(**row))
summary += [
    "",
    "C1 is partially supported: passed=%s; not passed=%s." % (passed, failed),
    "C2 direction-only descriptive replication (both expected negative): direction_consistent=%s; no registered model or significance threshold." % bool(c2_pass),
    "C3 comprises 20 items x 5 option positions, 4,800 responses; minimum cell n=2; descriptive only, with no registered sparsity failure threshold.",
    "The batch-1 versus batch-2 mean-scale table is post hoc and is not registered C1.",
    "E3a and order analyses are protocol-specified but unregistered and non-confirmatory.",
    "No item text, mapping, scoring key, parameter, or backend formula is disclosed.",
]
with open(os.path.join(OUT, "E3b_分析摘要.md"), "w", encoding="utf-8") as handle:
    handle.write("\n".join(summary) + "\n")
print("C1 passed", passed, "failed", failed, "C2 directions passed", bool(c2_pass))
