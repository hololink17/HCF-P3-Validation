# a24_item_equivalence.py —— 公开读数分布、外部题目统计与症状模式
# 内部题目级分支仅在受控输入存在时运行；公开 L2 不提供这些输入。
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import common as C

MODULE = "a24_item_equivalence"
ITEMS = ["hcf_%d" % i for i in range(5, 25)]


def load_items(d):
    cols = [c for c in ITEMS if c in d.columns]
    return d[cols].dropna()


def cross_cohort_equivalence(p, j):
    rows = []
    for c in ITEMS:
        if c not in p.columns or c not in j.columns:
            continue
        a = p[c].dropna(); b = j[c].dropna()
        if len(a) < 3 or len(b) < 3:
            continue
        t, pv = ttest_ind(a, b, equal_var=False)
        rows.append((c, round(a.mean(), 3), round(b.mean(), 3),
                     round(a.mean() - b.mean(), 3), round(t, 4), round(pv, 5), len(a), len(b)))
    out = pd.DataFrame(rows, columns=["题目", "Prolific均值", "见数均值", "差", "t", "p", "n_P", "n_J"])
    return out


def symptom_pattern(d, cohort):
    # 按 PHQ-9 总分区间描述 H 分布；不是同分匹配、条件独立或等价性检验。
    sub = d[["phq9_total", "H_score"]].dropna()
    # 以 phq9_total 分层（0-4 / 5-9 / 10-14 / 15+）
    bins = [-0.1, 4.5, 9.5, 14.5, 99]
    labs = ["0-4", "5-9", "10-14", "15+"]
    sub = sub.copy()
    sub["phq_band"] = pd.cut(sub["phq9_total"], bins=bins, labels=labs)
    grp = sub.groupby("phq_band", observed=True)["H_score"].agg(["mean", "std", "count"]).round(4)
    grp = grp.reset_index()
    grp.insert(0, "队列", cohort)
    return grp


def main():
    p = C.load_prolific(); j = C.load_jianshu()
    if any(c in p.columns and c in j.columns for c in ITEMS):
        C.save_tab(cross_cohort_equivalence(p, j), "a24_题目跨队列等价", "prolific", MODULE)
        C.save_tab(cross_cohort_equivalence(p, j), "a24_题目跨队列等价", "jianshu", MODULE)
    else:
        print("[a24] 跳过题目等价：公开 L2 不含受控题目级输入")
    for cohort, d in [("prolific", p), ("jianshu", j)]:
        C.save_tab(symptom_pattern(d, cohort), "a24_症状模式_H分层", cohort, MODULE)
    print("[a24] 公开症状模式完成")


if __name__ == "__main__":
    main()
