# a12_convergent.py —— 聚合效度 / MTMM / 外部验证 / 可接受性
# 移植 R 10_analysis_12_convergent_semantic.R + 22_incremental_mtmm.R + 20_external_validation.R + 21_acceptability.R
# 金标准收敛（n=142）：α/H/载体 与 MAIA/AAQ/DCS/BADS 的相关（聚合效度）。
import numpy as np
import pandas as pd
import common as C

MODULE = "a12_convergent"
GOLD_SCORES = {
    "maia_noticing": "MAIA-2 觉察", "maia_not_distracting": "MAIA-2 不分散(反向)",
    "maia_not_worrying": "MAIA-2 不担忧(部分反向)", "maia_three_dim": "MAIA-2 三维均",
    "aaq_total": "AAQ-II 经验性回避(越高越不灵活)", "dcs_total": "DCS 决策困难",
    "bads_total": "BADS 行为激活",
}
PRED = ["alpha1", "alpha2", "alpha3", "alpha4", "carrier_a", "carrier_b",
        "carrier_c", "carrier_d", "H_score"]


def convergent(d, cohort):
    rows = []
    for gcol, glabel in GOLD_SCORES.items():
        if gcol not in d.columns:
            continue
        for pred in PRED:
            if pred not in d.columns:
                continue
            sub = d[[gcol, pred]].dropna()
            if len(sub) < 10:
                continue
            r, p = C.pearson_r(sub[gcol], sub[pred])
            rows.append((glabel, gcol, pred, round(r, 4), round(p, 5), len(sub)))
    out = pd.DataFrame(rows, columns=["金标准", "金标准列", "HCF预测子", "r", "p", "n"])
    out.insert(0, "队列", cohort)
    return out


def acceptability(d, cohort):
    rows = []
    if "duration_sec" in d.columns:
        dur = d["duration_sec"].dropna()
        rows.append(("作答时长_均值秒", round(dur.mean(), 1)))
        rows.append(("作答时长_中位数秒", round(dur.median(), 1)))
    if "attention_pass" in d.columns:
        ap = d["attention_pass"].dropna()
        rows.append(("注意力通过率", round(ap.mean(), 4)))
    if "inattentive" in d.columns:
        ia = d["inattentive"].dropna()
        rows.append(("误答比例", round(ia.mean(), 4)))
    out = pd.DataFrame(rows, columns=["指标", "数值"])
    out.insert(0, "队列", cohort)
    return out


def main():
    g = C.load_gold()
    if g is not None:
        C.save_tab(convergent(g, "Prolific"), "a12_金标准收敛_n142", "prolific", MODULE)
        print("[a12] 金标准收敛（n=%d）完成" % len(g))
    else:
        print("[a12] 无金标准子样本")
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        acc = acceptability(d, cohort)
        if acc.empty:
            print("[a12] %s：公开 L2 不含作答时长或注意力检查字段；可接受性分析显式跳过，不生成占位结果" % cohort)
        else:
            C.save_tab(acc, "a12_可接受性", cohort, MODULE)
    print("[a12] 聚合效度完成；可接受性按输入字段处理")


if __name__ == "__main__":
    main()
