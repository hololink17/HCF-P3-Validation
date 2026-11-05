# a02_alpha_structure.py —— α 结构 / 七因子 PCA / 题目 PCA
# 移植 R 10_analysis_02_alpha_structure.R + 26_item_structure_efa.R + 30_nesting_structure.R
import numpy as np
import pandas as pd
import common as C

MODULE = "a02_alpha_structure"
SEVEN = ["alpha1", "alpha2", "alpha4", "carrier_a", "carrier_b", "carrier_c", "carrier_d"]


def seven_factor_pca(d, cohort):
    X = d[SEVEN].dropna().copy()
    Xz = (X - X.mean()) / X.std(ddof=0)
    R = np.corrcoef(Xz.values, rowvar=False)
    evals, evecs = np.linalg.eigh(R)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    var_ratio = evals / evals.sum()
    loadings = evecs[:, 0]  # PC1
    rows = []
    for i, name in enumerate(SEVEN):
        rows.append((name, round(loadings[i], 4),
                     round(np.corrcoef(Xz[name], Xz.dot(loadings))[0, 1], 4)))
    pc1 = pd.DataFrame(rows, columns=["因子", "PC1载荷", "与PC1相关"])
    pc1.insert(0, "队列", cohort)
    summ = pd.DataFrame([
        ("PC1方差解释比", round(var_ratio[0], 4)),
        ("PC2方差解释比", round(var_ratio[1], 4)),
        ("PC1+PC2", round(var_ratio[:2].sum(), 4)),
        ("七因子两两相关最小", round(np.min(R[np.triu_indices_from(R, 1)]), 4)),
        ("七因子两两相关最大", round(np.max(R[np.triu_indices_from(R, 1)]), 4)),
    ], columns=["指标", "数值"])
    summ.insert(0, "队列", cohort)
    return pc1, summ


def alpha_cronbach(d, cohort):
    items = d[["alpha1", "alpha2", "alpha3", "alpha4"]].dropna()
    a = C.cronbach_alpha(items.values)
    corr = items.corr().round(4)
    corr.insert(0, "因子", corr.columns)
    corr.insert(0, "队列", cohort)
    return pd.DataFrame([("α1-α4 Cronbach α", round(a, 4))], columns=["指标", "数值"]).assign(队列=cohort), corr


def pca_items(d, cohort):
    item_cols = ["hcf_%d" % i for i in range(5, 25) if "hcf_%d" % i in d.columns]
    if not item_cols:
        return None
    X = d[item_cols].dropna()
    Xz = (X - X.mean()) / X.std(ddof=0)
    R = np.corrcoef(Xz.values, rowvar=False)
    evals, evecs = np.linalg.eigh(R)
    order = np.argsort(evals)[::-1]
    evals = evals[order]
    out = pd.DataFrame({
        "成分": ["PC%d" % (i + 1) for i in range(len(evals))],
        "特征值": np.round(evals, 4),
        "方差解释比": np.round(evals / evals.sum(), 4),
        "累计": np.round(np.cumsum(evals) / evals.sum(), 4),
    })
    out.insert(0, "队列", cohort)
    return out


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        pc1, summ = seven_factor_pca(d, cohort)
        C.save_tab(pc1, "a02_七因子PC1", cohort, MODULE)
        C.save_tab(summ, "a02_七因子摘要", cohort, MODULE)
        a_row, corr = alpha_cronbach(d, cohort)
        C.save_tab(a_row, "a02_alpha_cronbach", cohort, MODULE)
        C.save_tab(corr, "a02_alpha相关阵", cohort, MODULE)
        item_pca = pca_items(d, cohort)
        if item_pca is None:
            print("[a02] %s：公开 L2 不含题目级作答；受控 PCA 显式跳过，不生成占位结果" % cohort)
        else:
            C.save_tab(item_pca, "a02_20题PCA", cohort, MODULE)
        # 分型分布
        prof = d["alpha_group"].value_counts(dropna=False).reset_index()
        prof.columns = ["alpha_group", "n"]
        prof["pct"] = (100 * prof["n"] / len(d)).round(1)
        prof.insert(0, "队列", cohort)
        C.save_tab(prof, "a02_分型分布", cohort, MODULE)
    print("[a02] α 结构 / 七因子完成；题目级 PCA 按输入层级处理")


if __name__ == "__main__":
    main()
