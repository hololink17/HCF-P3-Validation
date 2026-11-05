# a27_suicidal_ideation.py —— PHQ-9第9条（自杀意念）及α1增量效度
# 题项必须由构建层按官方0–3计分；缺失不得转成0。
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2
from scipy import stats
import common as C

MODULE = "a27_suicidal_ideation"


def _fit(d, y, preds):
    dd = d[[y] + preds].dropna()
    if len(dd) < 20 or dd[y].nunique() < 2:
        return None, dd
    X = sm.add_constant(dd[preds].astype(float), has_constant="add")
    try:
        fit = sm.GLM(dd[y].astype(float), X, family=sm.families.Binomial()).fit()
        return (fit if bool(fit.converged) else None), dd
    except Exception:
        return None, dd


def _incremental(d, base, added, label, cohort):
    m0, d0 = _fit(d, "si", base)
    m1, d1 = _fit(d, "si", base + added)
    if m0 is None or m1 is None or len(d0) != len(d1):
        return [cohort, label, len(d1), np.nan, np.nan, np.nan, np.nan]
    lr = max(0.0, 2 * (m1.llf - m0.llf))
    p = chi2.sf(lr, len(added))
    return [cohort, label, len(d1), m1.params[added[0]], np.exp(m1.params[added[0]]),
            m1.pvalues[added[0]], p]


def _item_profiles(d, cohort):
    """四α条目描述、组间检验及控制其余8题的偏相关；只输出聚合统计。"""
    rows = []
    other8 = d[["phq_%d" % i for i in range(36, 44)]].sum(axis=1, min_count=8)
    for alpha in ["alpha1", "alpha2", "alpha3", "alpha4"]:
        sub = pd.DataFrame({"item": d["phq_44"], "si": d["si"],
                            "alpha": d[alpha], "other8": other8}).dropna()
        pos = sub.loc[sub["si"] == 1, "alpha"].to_numpy(float)
        neg = sub.loc[sub["si"] == 0, "alpha"].to_numpy(float)
        pooled = np.sqrt(((len(pos)-1)*pos.var(ddof=1)+(len(neg)-1)*neg.var(ddof=1)) /
                         (len(pos)+len(neg)-2))
        effect = (pos.mean()-neg.mean())/pooled
        rho, rho_p = stats.spearmanr(sub["alpha"], sub["item"])
        x = sm.add_constant(sub[["other8"]].astype(float), has_constant="add")
        ra = sm.OLS(sub["alpha"].astype(float), x).fit().resid
        ri = sm.OLS(sub["item"].astype(float), x).fit().resid
        r = float(stats.pearsonr(ra, ri).statistic)
        df = len(sub)-3  # 一个控制变量：其余8题总分，不含第9题本身
        t = r*np.sqrt(df/(1-r*r))
        partial_p = float(2*stats.t.sf(abs(t), df))
        rows.append({"队列": cohort, "alpha": alpha, "n": len(sub),
                     "SI阳性数": int(sub["si"].sum()), "SI阳性率": sub["si"].mean(),
                     "阳性组均值": pos.mean(), "阴性组均值": neg.mean(),
                     "均值差_阳性减阴性": pos.mean()-neg.mean(), "Cohen_d": effect,
                     "Welch_p": stats.ttest_ind(pos, neg, equal_var=False).pvalue,
                     "MannWhitney_p": stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue,
                     "Spearman_r": rho, "Spearman_p": rho_p,
                     "partial_r_控制其余8题": r, "partial_df": df, "partial_p": partial_p})
    return pd.DataFrame(rows)


def _run(d, cohort):
    if "phq_44" not in d.columns:
        return
    d = d.copy()
    d["si"] = np.where(d["phq_44"].notna(), (d["phq_44"] >= 1).astype(float), np.nan)
    item = d[["phq_44", "si", "alpha1", "H_score", "phq9_total"]].dropna()
    C.save_tab(_item_profiles(d, cohort), "a27_条目描述与四alpha", cohort, MODULE)
    rows = []
    for pred in ["alpha1", "H_score", "phq9_total"]:
        if pred == "phq9_total":
            r, p = C.pearson_r(item[pred], item["si"])
        else:
            fit, dd = _fit(d, "si", [pred])
            r, p = (np.exp(fit.params[pred]), fit.pvalues[pred]) if fit is not None else (np.nan, np.nan)
        rows.append((pred, r, p, len(item)))
    r_s, p_s = C.spearman_r(item["alpha1"], item["phq_44"])
    rows.append(("alpha1_vs_item_spearman_r", r_s, p_s, len(item)))
    tab = pd.DataFrame(rows, columns=["预测子或指标", "OR或r", "p", "n"])
    tab.insert(0, "队列", cohort)
    C.save_tab(tab, "a27_自杀意念", cohort, MODULE)

    phq_other = ["phq_%d" % i for i in range(36, 44)]
    inc = []
    if all(c in d for c in phq_other):
        d["phq_other8"] = d[phq_other].sum(axis=1, min_count=8)
        inc.append(_incremental(d, ["phq_other8"], ["alpha1"], "PHQ其余8题+alpha1", cohort))
    if "maia_three_dim" in d and d["maia_three_dim"].notna().sum() >= 20:
        inc.append(_incremental(d, ["maia_three_dim"], ["alpha1"], "MAIA三维+alpha1", cohort))
    C.save_tab(pd.DataFrame(inc, columns=["队列", "模型", "n", "alpha1_logOR",
                                         "alpha1_OR", "alpha1_p", "增量LR_p"]),
               "a27_alpha1增量效度", cohort, MODULE)
    print("[a27] %s n=%d n(SI+)=%d Spearman(alpha1,item9)=%.4f" %
          (cohort, len(item), int(item["si"].sum()), r_s))


def main():
    for cohort in C.COHORTS:
        _run(C.load_cohort(cohort), cohort)
    gold = C.load_gold()
    if gold is not None:
        _run(gold, "gold")
    print("[a27] 自杀意念条目完成")


if __name__ == "__main__":
    main()
