# a19_ibs_subgroup.py —— IBS 亚组敏感性（移植 R 10_analysis_19_ibs_subgroup_sensitivity.R）
# 在性别 / 年龄亚组内，H_score 对 IBS 的判别力 AUC 是否稳定。
import numpy as np
import pandas as pd
import common as C

MODULE = "a19_ibs_subgroup"


def subgroup(d, cohort):
    rows = []
    y = d["ibs_positive"].dropna()
    s = d["H_score"].dropna()
    base = y.reindex(s.index).notna() & s.notna()
    # 性别
    for g in d["gender"].dropna().unique():
        m = base & (d["gender"] == g)
        if m.sum() < 20 or y[m].nunique() < 2:
            continue
        a, lo, hi = C.discrimination(y[m].values, s[m].values)
        rows.append(("性别=%s" % g, round(a, 4), round(lo, 4), round(hi, 4), int(m.sum())))
    # 年龄
    for alevel in d["age_group"].dropna().unique():
        m = base & (d["age_group"] == alevel)
        if m.sum() < 20 or y[m].nunique() < 2:
            continue
        a, lo, hi = C.discrimination(y[m].values, s[m].values)
        rows.append(("年龄=%s" % alevel, round(a, 4), round(lo, 4), round(hi, 4), int(m.sum())))
    out = pd.DataFrame(rows, columns=["亚组", "AUC", "CI低", "CI高", "n"])
    out.insert(0, "队列", cohort)
    return out


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        C.save_tab(subgroup(d, cohort), "a19_IBS亚组", cohort, MODULE)
        print("[a19] %s IBS 亚组敏感性完成" % cohort)
    print("[a19] IBS 亚组敏感性完成")


if __name__ == "__main__":
    main()
