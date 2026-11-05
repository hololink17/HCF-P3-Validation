# a01_descriptives.py —— 描述统计（移植 R 10_analysis_01_descriptives.R）
# 两队列：样本量、人口学分布、α/carrier/H 描述、结局患病率、金标准均值。
import numpy as np
import pandas as pd
import common as C

MODULE = "a01_descriptives"


def cohort_descriptives(d, cohort):
    rows = []
    n = len(d)
    rows.append(("样本量 n", n))
    # 人口学
    for col, label in [("age_group", "年龄组"), ("gender", "性别"), ("education", "教育")]:
        if col in d.columns:
            vc = d[col].value_counts(dropna=False)
            for k, v in vc.items():
                rows.append(("%s=%s" % (label, k), v, round(100 * v / n, 1)))
    # α/carrier/H 描述
    for col in ["alpha1", "alpha2", "alpha3", "alpha4",
                "carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score"]:
        if col in d.columns:
            s = d[col].dropna()
            rows.append((col + "_mean", round(s.mean(), 4)))
            rows.append((col + "_sd", round(s.std(ddof=1), 4)))
            rows.append((col + "_median", round(s.median(), 4)))
    # 结局患病率
    for col, label in [("depression", "抑郁筛查阳性(PHQ≥10)"),
                       ("ibs_positive", "IBS阳性(Rome IV)"),
                       ("chronic_pain", "慢性痛"), ("heart_disease", "心脏病"),
                       ("hypertension", "高血压"), ("diabetes", "糖尿病"),
                       ("hyperlipidemia", "高血脂"), ("obese", "肥胖(BMI≥30)")]:
        if col in d.columns:
            s = d[col].dropna()
            if s.nunique() <= 2:
                pos = int((s == 1).sum())
                rows.append((label + "_n", pos, round(100 * pos / len(s), 1)))
    out = pd.DataFrame(rows, columns=["指标", "数值", "百分比"]).fillna("")
    out.insert(0, "队列", cohort)
    return out


def gold_descriptives(g):
    rows = []
    for col in ["maia_noticing", "maia_not_distracting", "maia_not_worrying",
                "maia_three_dim", "aaq_total", "dcs_total", "bads_total"]:
        if col in g.columns:
            s = g[col].dropna()
            rows.append((col, round(s.mean(), 3), round(s.std(ddof=1), 3), len(s)))
    out = pd.DataFrame(rows, columns=["金标准量表", "均值", "标准差", "n"])
    return out


def main():
    results = []
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        results.append(cohort_descriptives(d, cohort))
    combined = pd.concat(results, ignore_index=True)
    C.save_tab(combined, "a01_描述统计_按队列", "prolific", MODULE)
    C.save_tab(combined, "a01_描述统计_按队列", "jianshu", MODULE)
    g = C.load_gold()
    if g is not None:
        gd = gold_descriptives(g)
        C.save_tab(gd, "a01_金标准描述_n142", "prolific", MODULE)
    print("[a01] 描述统计完成")
    print(combined.to_string(index=False))


if __name__ == "__main__":
    main()
