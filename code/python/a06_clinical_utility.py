# a06_clinical_utility.py —— 临床效用 / PPV-NPV（移植 R 10_analysis_06_clinical_utility.R + 23_clinical_operability.R）
# 以 H_score 低分→高风险，对抑郁/IBS 在多个 H 切点下计算敏感度/特异度/PPV/NPV/Youden。
import numpy as np
import pandas as pd
import common as C

MODULE = "a06_clinical_utility"


def utility(d, cohort, outcome, score="H_score", cuts=None):
    y = d[outcome].dropna()
    s = d[score].dropna()
    idx = y.index.intersection(s.index)
    y = y.loc[idx].astype(int).values
    s = s.loc[idx].values
    if cuts is None:
        cuts = np.nanquantile(s, [0.10, 0.20, 0.30, 0.40])
    base_rate = y.mean()
    rows = []
    for c in cuts:
        # H 低 → 风险高：预测阳性 = s <= c
        pred = (s <= c).astype(int)
        tp = int(((pred == 1) & (y == 1)).sum())
        fp = int(((pred == 1) & (y == 0)).sum())
        fn = int(((pred == 0) & (y == 1)).sum())
        tn = int(((pred == 0) & (y == 0)).sum())
        sens = tp / (tp + fn) if (tp + fn) else np.nan
        spec = tn / (tn + fp) if (tn + fp) else np.nan
        ppv = tp / (tp + fp) if (tp + fp) else np.nan
        npv = tn / (tn + fn) if (tn + fn) else np.nan
        youden = (sens if sens == sens else 0) + (spec if spec == spec else 0) - 1
        rows.append((round(float(c), 4), int((s <= c).sum()), tp, fp, fn, tn,
                     round(sens, 4), round(spec, 4), round(ppv, 4), round(npv, 4),
                     round(youden, 4)))
    out = pd.DataFrame(rows, columns=["H切点", "预测阳性数", "TP", "FP", "FN", "TN",
                                      "敏感度", "特异度", "PPV", "NPV", "Youden"])
    out.insert(0, "结局", outcome)
    out.insert(0, "队列", cohort)
    return out, base_rate


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        for outcome in ["depression", "ibs_positive"]:
            if outcome not in d.columns:
                continue
            out, br = utility(d, cohort, outcome)
            C.save_tab(out, "a06_临床效用_%s" % outcome, cohort, MODULE)
            print("[a06] %s %s 基础患病率=%.3f" % (cohort, outcome, br))
    print("[a06] 临床效用/PPV-NPV 完成")


if __name__ == "__main__":
    main()
