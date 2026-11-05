# a18_anchor_vs_algorithm.py —— 锚定题 vs 算法（移植 R 10_analysis_18_anchor_vs_algorithm.R）
# 锚定题（A=5…F=0，四题均 A=最健康，无反向）与算法 α1-α4 的相关 + 一致性。
import numpy as np
import pandas as pd
import common as C

MODULE = "a18_anchor_vs_algorithm"
ANCHOR = ["anchor1", "anchor2", "anchor3", "anchor4"]
ALPHA = ["alpha1", "alpha2", "alpha3", "alpha4"]


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        rows = []
        for a, al in zip(ANCHOR, ALPHA):
            sub = d[[a, al]].dropna()
            if len(sub) < 10:
                continue
            r, p = C.pearson_r(sub[a], sub[al])
            rows.append((a, al, round(r, 4), round(p, 5), len(sub)))
        tab = pd.DataFrame(rows, columns=["锚定题", "算法因子", "Pearson_r", "p", "n"])
        tab.insert(0, "队列", cohort)
        C.save_tab(tab, "a18_锚定vs算法相关", cohort, MODULE)
        # 锚定题内部一致性
        anc = d[ANCHOR].dropna()
        a_alpha = C.cronbach_alpha(anc.values)
        C.save_tab(pd.DataFrame([("锚定四题 Cronbach α", round(a_alpha, 4)),
                                 ("锚定均值", round(anc.mean().mean(), 4))],
                               columns=["指标", "数值"]).assign(队列=cohort),
                  "a18_锚定可靠性", cohort, MODULE)
        print("[a18] %s 锚定题 vs 算法完成 | α=%.3f" % (cohort, a_alpha))
    print("[a18] 锚定题 vs 算法完成")


if __name__ == "__main__":
    main()
