# a15_screening_roc.py —— 筛查 ROC 与最优切点（移植 R 10_analysis_15_screening_roc.R）
# 对抑郁/IBS，以 H_score（及 α1/α4）为预测子，ROC 最优切点（Youden）、敏感度/特异度/PPV/NPV/AUC。
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve
import common as C

MODULE = "a15_screening_roc"
PRED = ["H_score", "alpha1", "alpha4"]


def roc_cut(d, cohort, outcome):
    y = d[outcome].dropna()
    rows = []
    for pred in PRED:
        s = d[pred].dropna()
        m = y.reindex(s.index).notna() & s.notna()
        if m.sum() < 20 or y[m].nunique() < 2:
            continue
        yy = y[m].values.astype(int)
        ss = s[m].values
        auc, lo, hi = C.discrimination(yy, ss)
        # Youden 最优切点（H 越低风险越高 → 以 -ss 做 ROC）
        fpr, tpr, thr = roc_curve(yy, -ss)
        youden = tpr - fpr
        k = int(np.argmax(youden))
        cut = -thr[k]  # 还原为 H 的原刻度切点
        pred_pos = (ss <= cut).astype(int)
        tp = int(((pred_pos == 1) & (yy == 1)).sum()); fp = int(((pred_pos == 1) & (yy == 0)).sum())
        fn = int(((pred_pos == 0) & (yy == 1)).sum()); tn = int(((pred_pos == 0) & (yy == 0)).sum())
        sens = tp / (tp + fn) if (tp + fn) else np.nan
        spec = tn / (tn + fp) if (tn + fp) else np.nan
        ppv = tp / (tp + fp) if (tp + fp) else np.nan
        npv = tn / (tn + fn) if (tn + fn) else np.nan
        rows.append((outcome, pred, round(auc, 4), round(lo, 4), round(hi, 4),
                     round(float(cut), 4), round(sens, 4), round(spec, 4),
                     round(ppv, 4), round(npv, 4)))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "AUC", "CI低", "CI高",
                                      "Youden切点", "敏感度", "特异度", "PPV", "NPV"])
    out.insert(0, "队列", cohort)
    return out


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        for outcome in ["depression", "ibs_positive"]:
            if outcome not in d.columns:
                continue
            C.save_tab(roc_cut(d, cohort, outcome), "a15_ROC筛查_%s" % outcome, cohort, MODULE)
        print("[a15] %s ROC 筛查完成" % cohort)
    print("[a15] 筛查 ROC 与最优切点完成")


if __name__ == "__main__":
    main()
