# a10_robustness.py —— 稳健性 / 敏感性（移植 R 10_analysis_10_robustness.R + 13_sensitivity.R）
# 判别力 AUC 在年龄/性别分层下的稳定性；剔除注意力失败者后的 AUC 变化。
# 公开 L2 不含注意力字段；教育分层与 bootstrap CI 不在本模块的实际输出中。
import numpy as np
import pandas as pd
import common as C

MODULE = "a10_robustness"
PRED = ["H_score", "alpha1", "alpha4"]
OUTC = ["depression", "ibs_positive"]


def auc_strata(d, cohort):
    rows = []
    for outcome in OUTC:
        if outcome not in d.columns:
            continue
        y = d[outcome].dropna()
        for pred in PRED:
            s = d[pred].dropna()
            m = y.reindex(s.index).notna() & s.notna()
            base_auc, _, _ = C.discrimination(y[m].values, s[m].values)
            rows.append((outcome, pred, "全体", round(base_auc, 4), int(m.sum())))
            # 性别分层
            for gval in d["gender"].dropna().unique():
                gi = (d["gender"] == gval) & m
                if gi.sum() < 20:
                    continue
                a, _, _ = C.discrimination(y[gi].values, s[gi].values)
                rows.append((outcome, pred, "性别=%s" % gval, round(a, 4), int(gi.sum())))
            # 年龄分层
            for alevel in d["age_group"].dropna().unique():
                ai = (d["age_group"] == alevel) & m
                if ai.sum() < 20:
                    continue
                a2, _, _ = C.discrimination(y[ai].values, s[ai].values)
                rows.append((outcome, pred, "年龄=%s" % alevel, round(a2, 4), int(ai.sum())))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "分层", "AUC", "n"])
    out.insert(0, "队列", cohort)
    C.save_tab(out, "a10_分层稳健性", cohort, MODULE)


def attention_sensitivity(d, cohort):
    rows = []
    if "inattentive" not in d.columns and "attention_pass" not in d.columns:
        print("[a10] %s：公开 L2 不含注意力字段；注意力敏感性显式跳过，不生成占位结果" % cohort)
        return
    for outcome in OUTC:
        if outcome not in d.columns:
            continue
        flag = d["inattentive"].fillna(0) if "inattentive" in d.columns else (1 - d["attention_pass"].fillna(1))
        for pred in PRED:
            s = d[pred].dropna()
            y = d.loc[s.index, outcome].dropna()
            m = y.notna() & s.notna()
            a_full, _, _ = C.discrimination(y[m].values, s[m].values)
            keep = m & (flag == 0)
            a_clean, _, _ = C.discrimination(y[keep].values, s[keep].values) if keep.sum() > 20 else (np.nan, np.nan, np.nan)
            rows.append((outcome, pred, round(a_full, 4), round(a_clean, 4),
                         int(m.sum()), int(keep.sum())))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "AUC_全体", "AUC_剔除误答", "n_全体", "n_清洁"])
    out.insert(0, "队列", cohort)
    C.save_tab(out, "a10_注意力敏感性", cohort, MODULE)


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        auc_strata(d, cohort)
        attention_sensitivity(d, cohort)
        print("[a10] %s 分层分析完成；注意力分支按输入层级处理" % cohort)
    print("[a10] 分层分析完成；注意力分支未执行时不得计为敏感性证据")


if __name__ == "__main__":
    main()
