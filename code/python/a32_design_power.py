# a32_design_power.py —— 设计层样本量与功效（移植 R 10_analysis_14_power.R + 32_design_power.R）
# ① 已实现的判别力功效（AUC≠0.5 的 post-hoc power）；② 预注册设计参数文档
#    （E3a 重测 n=160 证 ICC；E3b TOST 等效 |δ|<0.147，d=0.262 为界）。
import numpy as np
import pandas as pd
from scipy.stats import norm
import common as C

MODULE = "a32_design_power"
PRED = ["H_score", "alpha1", "alpha4"]
OUTC = ["depression", "ibs_positive"]


def posthoc_power(auc, se):
    if se <= 0 or not np.isfinite(auc):
        return np.nan
    z = (auc - 0.5) / se
    # 与设计表声明的双侧 alpha=0.05 一致；仍仅为观测功效描述，不替代前瞻设计。
    return float(norm.cdf(z - 1.96) + norm.cdf(-z - 1.96))


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        rows = []
        for outcome in OUTC:
            if outcome not in d.columns:
                continue
            y = d[outcome].dropna()
            for pred in PRED:
                s = d[pred].dropna()
                m = y.reindex(s.index).notna() & s.notna()
                if m.sum() < 20 or y[m].nunique() < 2:
                    continue
                a, lo, hi = C.discrimination(y[m].values, s[m].values)
                se = (hi - lo) / (2 * 1.96)
                pw = posthoc_power(a, se)
                rows.append((outcome, pred, round(a, 4), round(se, 4),
                             round((hi - lo) / 2, 4), round(pw, 4), int(m.sum())))
        tab = pd.DataFrame(rows, columns=["结局", "预测因子", "AUC", "SE", "CI半宽", "post-hoc功效", "n"])
        tab.insert(0, "队列", cohort)
        C.save_tab(tab, "a32_功效", cohort, MODULE)
        print("[a32] %s 设计功效完成" % cohort)
    # 预注册设计参数文档
    design = pd.DataFrame([
        ("E3a 重测信度 ICC", "目标 n=160（同批人同题目不同时间）", "证 ICC≥0.7"),
        ("E3b 跨现象域独立同构", "TOST 等效 |δ|<0.147 (d=0.262)", "四因子结构一致性"),
        ("本研究实际 n", "Prolific 462 / 见数 504", "高于最小样本要求"),
        ("判别力检验 α", "0.05 双侧", "—"),
    ], columns=["设计项", "参数", "说明"])
    C.save_tab(design, "a32_预注册设计参数", "prolific", MODULE)
    C.save_tab(design, "a32_预注册设计参数", "jianshu", MODULE)
    print("[a32] 设计层样本量与功效完成")


if __name__ == "__main__":
    main()
