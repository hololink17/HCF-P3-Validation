# a11_cross_cohort.py —— 跨队列对照（移植 R 10_analysis_11_cross_cohort.R + 20_external_validation.R）
# Prolific vs 见数：α/H/载体/结局均值差异（Welch t）；抑郁/IBS 判别力 AUC 差异 bootstrap CI。
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import common as C

MODULE = "a11_cross_cohort"
COLS = ["alpha1", "alpha2", "alpha3", "alpha4", "carrier_a", "carrier_b",
        "carrier_c", "carrier_d", "H_score", "phq9_total", "gad7_total",
        "ibs_symptom_count", "ibs_positive", "depression"]


def main():
    p = C.load_prolific()
    j = C.load_jianshu()
    # 均值对照
    rows = []
    for c in COLS:
        a = p[c].dropna(); b = j[c].dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        t, pv = ttest_ind(a, b, equal_var=False)
        rows.append((c, round(a.mean(), 4), round(b.mean(), 4),
                     round(a.mean() - b.mean(), 4), round(t, 4), round(pv, 6)))
    comp = pd.DataFrame(rows, columns=["变量", "Prolific均值", "见数均值", "差", "t", "p"])
    C.save_tab(comp, "a11_均值对照", "prolific", MODULE)
    C.save_tab(comp, "a11_均值对照", "jianshu", MODULE)

    # AUC 差异 bootstrap（抑郁 / IBS）
    rows2 = []
    for outcome in ["depression", "ibs_positive"]:
        y1 = p[outcome].dropna(); s1 = p["H_score"].dropna()
        y2 = j[outcome].dropna(); s2 = j["H_score"].dropna()
        m1 = y1.reindex(s1.index).notna() & s1.notna()
        m2 = y2.reindex(s2.index).notna() & s2.notna()
        a1, _, _ = C.discrimination(y1[m1].values, s1[m1].values)
        a2, _, _ = C.discrimination(y2[m2].values, s2[m2].values)
        dlt, lo, hi, nb = C.bootstrap_diff(y1[m1].values, s1[m1].values,
                                           y2[m2].values, s2[m2].values)
        rows2.append((outcome, round(a1, 4), round(a2, 4), round(dlt, 4),
                      round(lo, 4), round(hi, 4), int(nb)))
    aucd = pd.DataFrame(rows2, columns=["结局", "Prolific_AUC", "见数_AUC",
                                        "AUC差", "CI低", "CI高", "bootstrap_n"])
    C.save_tab(aucd, "a11_AUC跨队列差异", "prolific", MODULE)
    C.save_tab(aucd, "a11_AUC跨队列差异", "jianshu", MODULE)
    print("[a11] 跨队列对照完成")
    print(comp.to_string(index=False))
    print(aucd.to_string(index=False))


if __name__ == "__main__":
    main()
