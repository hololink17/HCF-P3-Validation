# a05_ibs.py —— IBS 剂量反应 / 多切点（移植 R 10_analysis_05_ibs_dose_response.R + 19_ibs_subgroup_sensitivity.R + 35_cross_cohort_outcomes.R）
# IBS 症状计数 0-4 的剂量反应；H/α 在各症状数分组的均值；IBS阳性 vs 阴性 Cohen's d；
# 各切点（>=1/>=2/>=3/==4）判别力 AUC。
import numpy as np
import pandas as pd
from scipy.stats import f_oneway
import common as C

MODULE = "a05_ibs"
PRED = ["H_score", "alpha1", "alpha2", "alpha3", "alpha4"]


def dose_response(d, cohort):
    rows = []
    sc = d["ibs_symptom_count"].dropna()
    for k in range(5):
        m = sc == k
        if m.sum() == 0:
            continue
        row = {"队列": cohort, "症状数": k, "n": int(m.sum())}
        for p in PRED:
            s = d.loc[m, p].dropna()
            row[p] = round(s.mean(), 4) if len(s) else np.nan
        rows.append(row)
    dose = pd.DataFrame(rows)
    C.save_tab(dose, "a05_IBS剂量反应", cohort, MODULE)

    # ANOVA: H_score / α 跨症状数分组
    an = []
    for p in PRED:
        groups = [d.loc[d["ibs_symptom_count"] == k, p].dropna() for k in range(5)]
        groups = [g for g in groups if len(g) > 1]
        if len(groups) >= 2:
            f, pval = f_oneway(*groups)
            an.append((p, round(f, 4), round(pval, 6), len(groups)))
    anova = pd.DataFrame(an, columns=["预测子", "ANOVA_F", "p", "分组数"])
    anova.insert(0, "队列", cohort)
    C.save_tab(anova, "a05_IBS症状数_ANOVA", cohort, MODULE)


def multi_cut(d, cohort):
    rows = []
    sc = d["ibs_symptom_count"].dropna()
    for pred in PRED:
        s = d[pred].dropna()
        for cut, lab in [(1, ">=1"), (2, ">=2"), (3, ">=3"), (4, "==4")]:
            if cut == 4:
                y = (sc == 4).astype(int)
            else:
                y = (sc >= cut).astype(int)
                y = y.reindex(s.index)
            m = y.notna() & s.notna()
            if m.sum() < 20:
                continue
            auc, lo, hi = C.discrimination(y[m].values, s[m].values)
            rows.append((pred, lab, round(auc, 4), round(lo, 4), round(hi, 4),
                        int(m.sum()), int(y[m].sum())))
    out = pd.DataFrame(rows, columns=["预测因子", "IBS切点", "AUC", "CI低", "CI高", "n", "阳性数"])
    out.insert(0, "队列", cohort)
    C.save_tab(out, "a05_IBS多切点AUC", cohort, MODULE)


def ibs_vs_not(d, cohort):
    pos = d[d["ibs_positive"] == 1]
    neg = d[d["ibs_positive"] == 0]
    rows = []
    for p in PRED:
        a = pos[p].dropna().values; b = neg[p].dropna().values
        if len(a) > 3 and len(b) > 3:
            dd = C.cohens_d(np.concatenate([a, b]), np.r_[np.ones(len(a)), np.zeros(len(b))])
            rows.append((p, round(a.mean(), 4), round(b.mean(), 4), round(dd, 4)))
    out = pd.DataFrame(rows, columns=["预测因子", "IBS+均值", "IBS-均值", "Cohen_d(IBS+vs-)"])
    out.insert(0, "队列", cohort)
    C.save_tab(out, "a05_IBS阳性vs阴性", cohort, MODULE)


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        dose_response(d, cohort)
        multi_cut(d, cohort)
        ibs_vs_not(d, cohort)
        print("[a05] %s IBS 分析完成" % cohort)
    print("[a05] IBS 剂量反应/多切点完成")


if __name__ == "__main__":
    main()
