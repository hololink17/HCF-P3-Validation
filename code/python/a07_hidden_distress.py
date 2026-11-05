# a07_hidden_distress.py —— 潜在痛苦者分析（规范 A4，保留原定义并扩展外部效标）
# 主分析三组（互斥）：
#   健康对照 = 高意识型 & PHQ-9 < 10
#   潜在痛苦 = 耗竭型 & PHQ-9 < 10（不得额外排除 IBS/慢性痛等共病）
#   明确痛苦 = PHQ-9 >= 10
# 原“无任何临床阳性”的严格三组仅作敏感性分析，绝不用于检验慢性痛/IBS率差。
# 后台 H 与 α 分型同源，相关组间比较只作描述，不作独立外部效度证据。
import numpy as np
import pandas as pd
from scipy.stats import f_oneway, fisher_exact, ttest_ind
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import common as C

MODULE = "a07_hidden_distress"
CLINICAL = ["depression", "ibs_positive", "chronic_pain", "heart_disease",
            "hypertension", "diabetes", "hyperlipidemia"]
CRITERIA = ["phq9_total", "gad7_total", "ibs_symptom_count", "H_score",
            "alpha1", "alpha2", "alpha3", "alpha4"]
DISEASES = ["ibs_positive", "chronic_pain", "heart_disease", "hypertension",
            "diabetes", "hyperlipidemia", "obese", "sleep_disorder_positive"]
CONTINUOUS = ["phq_39", "sleep_hour", "sleep_symptom_count", "bmi", "H_score"]


def _cohens_d(x, g, a, b):
    xa = np.asarray(x)[np.asarray(g) == a]
    xb = np.asarray(x)[np.asarray(g) == b]
    xa = xa[np.isfinite(xa)]
    xb = xb[np.isfinite(xb)]
    if len(xa) < 3 or len(xb) < 3:
        return np.nan
    sp = np.sqrt(((len(xa) - 1) * xa.var(ddof=1) +
                  (len(xb) - 1) * xb.var(ddof=1)) / (len(xa) + len(xb) - 2))
    return (xa.mean() - xb.mean()) / sp if sp > 0 else np.nan


def assign_groups(d):
    """原始、不可替换的潜在痛苦主定义。"""
    g = np.full(len(d), "其他", dtype=object)
    hc = d["healthy_control"].fillna(False).values
    pd_ = d["potential_depression"].fillna(False).values
    manifest = (d["phq9_total"] >= 10).fillna(False).values
    g[manifest] = "明确痛苦"
    g[pd_ & ~manifest] = "潜在痛苦"
    g[hc & ~manifest] = "健康对照"
    return g


def assign_doubly_normal_groups(d):
    """方案口径敏感性：PHQ-9<10 且 GAD-7<10；不额外排除躯体共病。"""
    g = np.full(len(d), "其他", dtype=object)
    complete = d["phq9_total"].notna() & d["gad7_total"].notna()
    doubly_normal = complete & d["phq9_total"].lt(10) & d["gad7_total"].lt(10)
    manifest = d["phq9_total"].ge(10) | d["gad7_total"].ge(10)
    g[(d["alpha_group"] == "耗竭型") & doubly_normal] = "潜在痛苦"
    g[(d["alpha_group"] == "高意识型") & doubly_normal] = "健康对照"
    g[manifest.fillna(False)] = "明确痛苦"
    return g


def assign_strict_groups(d):
    """旧 Python 严格三组，仅保留作敏感性分析。"""
    any_clin = np.zeros(len(d), dtype=bool)
    clinical_complete = np.ones(len(d), dtype=bool)
    for c in CLINICAL:
        if c in d.columns:
            any_clin = any_clin | (d[c] == 1).fillna(False).values
            clinical_complete = clinical_complete & d[c].notna().values
    g = np.full(len(d), "其他", dtype=object)
    hc = d["healthy_control"].fillna(False).values
    pd_ = d["potential_depression"].fillna(False).values
    no_known_positive = ~any_clin & clinical_complete
    g[hc & no_known_positive] = "健康对照"
    g[pd_ & no_known_positive & ~hc] = "潜在痛苦"
    g[any_clin] = "明确痛苦"
    return g


def _size_table(d, group_col, cohort):
    size = d[group_col].value_counts().reindex(
        ["健康对照", "潜在痛苦", "明确痛苦", "其他"]).fillna(0).astype(int)
    out = size.reset_index()
    out.columns = ["组别", "n"]
    out["pct"] = (100 * out["n"] / len(d)).round(1)
    out.insert(0, "队列", cohort)
    return out


def _gradient(d, group_col, cohort):
    sub = d[d[group_col].isin(["健康对照", "潜在痛苦", "明确痛苦"])]
    gmap = {"健康对照": 0, "潜在痛苦": 1, "明确痛苦": 2}
    rows, tukey_rows = [], []
    for crit in CRITERIA:
        if crit not in sub.columns:
            continue
        s = sub[crit].dropna()
        grp = sub.loc[s.index, group_col].map(gmap)
        groups = [s[grp == k] for k in [0, 1, 2]]
        if min(map(len, groups)) < 2:
            continue
        f, p = f_oneway(*groups)
        tukey = pairwise_tukeyhsd(s.values, grp.values)
        for line in str(tukey).splitlines()[1:]:
            parts = line.split()
            if len(parts) >= 7 and parts[0] in {"0", "1", "2"}:
                tukey_rows.append([cohort, crit] + parts)
        rows.append((crit, round(groups[0].mean(), 4), round(groups[1].mean(), 4),
                     round(groups[2].mean(), 4), round(f, 4), round(p, 6),
                     round(_cohens_d(s.values, grp.values, 1, 0), 4),
                     round(_cohens_d(s.values, grp.values, 2, 0), 4)))
    grad = pd.DataFrame(rows, columns=["效标", "健康对照均值", "潜在痛苦均值",
                                      "明确痛苦均值", "ANOVA_F", "ANOVA_p",
                                      "d_潜在v健康", "d_明确v健康"])
    grad.insert(0, "队列", cohort)
    tk = pd.DataFrame(tukey_rows, columns=["队列", "效标", "组1", "组2", "均值差",
                                           "p_adj", "CI低", "CI高", "拒绝"])
    return sub, grad, tk


def _binary_comparison(a, b):
    a = pd.Series(a).dropna().astype(int)
    b = pd.Series(b).dropna().astype(int)
    tab = np.array([[int((a == 1).sum()), int((a == 0).sum())],
                    [int((b == 1).sum()), int((b == 0).sum())]])
    if (tab.sum(axis=1) == 0).any():
        return np.nan
    return fisher_exact(tab)[1]


def _disease_fatigue_table(d, cohort):
    """潜在痛苦核心价值：PHQ正常者中，耗竭型与高意识型的躯体/疾病差异。"""
    sub = d[d["group3"].isin(["健康对照", "潜在痛苦"])]
    hd = sub[sub["group3"] == "潜在痛苦"]
    hc = sub[sub["group3"] == "健康对照"]
    rows = []
    for v in DISEASES:
        if v not in sub.columns:
            continue
        a, b = hd[v].dropna(), hc[v].dropna()
        if len(a) == 0 or len(b) == 0:
            continue
        rows.append([cohort, v, "患病率", len(a), len(b), a.mean(), b.mean(),
                     a.mean() - b.mean(), _binary_comparison(a, b), np.nan])
    for v in CONTINUOUS:
        if v not in sub.columns:
            continue
        a, b = hd[v].dropna(), hc[v].dropna()
        if len(a) < 2 or len(b) < 2:
            continue
        p = ttest_ind(a, b, equal_var=False, nan_policy="omit").pvalue
        effect = _cohens_d(pd.concat([a, b]).values,
                           np.r_[np.ones(len(a)), np.zeros(len(b))], 1, 0)
        rows.append([cohort, v, "均值", len(a), len(b), a.mean(), b.mean(),
                     a.mean() - b.mean(), p, effect])
    return pd.DataFrame(rows, columns=["队列", "指标", "类型", "n_潜在痛苦", "n_健康对照",
                                      "潜在痛苦_均值或率", "健康对照_均值或率", "差值",
                                      "p_Fisher或Welch", "Cohen_d"])


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort).copy()
        d["group3"] = assign_groups(d)
        d["group3_strict_sensitivity"] = assign_strict_groups(d)

        C.save_tab(_size_table(d, "group3", cohort), "a07_三组规模", cohort, MODULE)
        sub, grad, tukey = _gradient(d, "group3", cohort)
        C.save_tab(grad, "a07_外部效标梯度", cohort, MODULE)
        C.save_tab(tukey, "a07_Tukey_HSD", cohort, MODULE)
        long = sub.melt(id_vars=["group3"], value_vars=[c for c in CRITERIA if c in sub],
                        var_name="效标", value_name="值").dropna()
        C.save_tab(long, "a07_效标长表", cohort, MODULE)
        C.save_tab(_disease_fatigue_table(d, cohort), "a07_潜在痛苦_疾病疲劳睡眠比较", cohort, MODULE)

        # 方案“双正常”口径与PHQ阈值敏感性均允许共病，只改变症状正常定义。
        d["group3_doubly_normal"] = assign_doubly_normal_groups(d)
        C.save_tab(_size_table(d, "group3_doubly_normal", cohort),
                   "a07_双正常三组规模_敏感性", cohort, MODULE)
        C.save_tab(_disease_fatigue_table(d.rename(columns={"group3": "group3_main"}).rename(
            columns={"group3_doubly_normal": "group3"}), cohort),
            "a07_双正常_疾病疲劳睡眠比较_敏感性", cohort, MODULE)
        threshold_tables = []
        for cutoff in [4, 6, 8, 10]:
            normal = d["phq9_total"].lt(cutoff)
            d["group3"] = "其他"
            d.loc[(d["alpha_group"] == "耗竭型") & normal, "group3"] = "潜在痛苦"
            d.loc[(d["alpha_group"] == "高意识型") & normal, "group3"] = "健康对照"
            tab = _disease_fatigue_table(d, cohort)
            tab.insert(1, "PHQ阈值", cutoff)
            threshold_tables.append(tab)
        C.save_tab(pd.concat(threshold_tables, ignore_index=True),
                   "a07_PHQ阈值4_6_8_10_敏感性", cohort, MODULE)
        d["group3"] = assign_groups(d)  # 阈值循环后恢复主三组，避免日志与后续调用误用末档分组

        # 严格定义只保留为敏感性分析；禁止从中解释疾病患病率。
        C.save_tab(_size_table(d, "group3_strict_sensitivity", cohort),
                   "a07_严格三组规模_仅敏感性", cohort, MODULE)
        _, strict_grad, _ = _gradient(d, "group3_strict_sensitivity", cohort)
        C.save_tab(strict_grad, "a07_严格三组效标梯度_仅敏感性", cohort, MODULE)
        note = pd.DataFrame([
            ["主定义", "潜在痛苦=耗竭型且PHQ-9<10；健康对照=高意识型且PHQ-9<10；允许并检验共病。"],
            ["严格组限制", "仅作敏感性分析；因额外排除全部临床阳性，不得用于疾病率差或慢性痛主张。"],
            ["H解释边界", "后台H与α分型同源，组间H差异仅描述，不是独立外部效度证据。"]
        ], columns=["项目", "说明"])
        C.save_tab(note, "a07_分析口径与边界", cohort, MODULE)
        print("[a07] %s 潜在痛苦主分析完成：\n%s" %
              (cohort, _size_table(d, "group3", cohort).to_string(index=False)))
    print("[a07] 潜在痛苦分析完成")


if __name__ == "__main__":
    main()
