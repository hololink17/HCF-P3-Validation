# a09_carrier.py —— 四载体系统分析（移植 R 10_analysis_09_carrier.R）
# 四个载体 a/b/c/d 的描述、相关、各自判别力 AUC、与 α 的关系。
import numpy as np
import pandas as pd
import common as C

MODULE = "a09_carrier"
CARRIERS = ["carrier_a", "carrier_b", "carrier_c", "carrier_d"]


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        # 描述
        desc = d[CARRIERS].describe().T.round(4).reset_index()
        desc.insert(0, "队列", cohort)
        C.save_tab(desc, "a09_载体描述", cohort, MODULE)
        # 相关
        R = d[CARRIERS].corr().round(4)
        C.save_tab(R, "a09_载体相关阵", cohort, MODULE)
        # 与 α 的相关
        Ra = d[CARRIERS + ["alpha1", "alpha2", "alpha3", "alpha4"]].corr().round(4)
        C.save_tab(Ra, "a09_载体与α相关", cohort, MODULE)
        # 判别力 AUC
        rows = []
        for outcome in ["depression", "ibs_positive", "chronic_pain", "hypertension"]:
            if outcome not in d.columns:
                continue
            y = d[outcome].dropna()
            for c in CARRIERS:
                s = d[c].dropna()
                m = y.reindex(s.index).notna() & s.notna()
                if m.sum() < 20:
                    continue
                auc, lo, hi = C.discrimination(y[m].values, s[m].values)
                rows.append((outcome, c, round(auc, 4), round(lo, 4), round(hi, 4)))
        tab = pd.DataFrame(rows, columns=["结局", "载体", "AUC", "CI低", "CI高"])
        tab.insert(0, "队列", cohort)
        C.save_tab(tab, "a09_载体判别力AUC", cohort, MODULE)
        print("[a09] %s 四载体分析完成" % cohort)
    print("[a09] 四载体系统分析完成")


if __name__ == "__main__":
    main()
