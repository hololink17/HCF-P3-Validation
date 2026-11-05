# a41_gold_incremental_ipsative.py —— Gold增量效度、症状偏相关与组内剖面
# 所有输入来自官方Gold基底(n=142)；只输出聚合统计，不输出逐人数据。
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

import common as C

MODULE = "a41_gold_incremental_ipsative"
ALPHAS = ["alpha1", "alpha2", "alpha3", "alpha4"]
GOLDS = ["maia_three_dim", "aaq_total", "dcs_total", "bads_total"]
SYMPTOMS = ["phq9_total", "gad7_total"]


def _design(d, cols):
    return sm.add_constant(d[cols].astype(float), has_constant="add")


def _ols_compare(d, outcome, base_cols, full_cols):
    sub = d[[outcome] + list(dict.fromkeys(base_cols + full_cols))].dropna()
    y = sub[outcome].astype(float)
    m0 = sm.OLS(y, _design(sub, base_cols)).fit()
    m1 = sm.OLS(y, _design(sub, full_cols)).fit()
    f, p, df_diff = m1.compare_f_test(m0)
    return {
        "outcome": outcome, "n": int(len(sub)), "R2_base": float(m0.rsquared),
        "R2_full": float(m1.rsquared), "delta_R2": float(m1.rsquared - m0.rsquared),
        "F": float(f), "df_numerator": float(df_diff), "df_denominator": float(m1.df_resid),
        "p": float(p), "status": "ok",
    }


def _paired_delong(y, score0, score1):
    """两相关ROC的DeLong检验，返回AUC0/AUC1/差值/z/p。"""
    y = np.asarray(y, dtype=int)
    scores = [np.asarray(score0, dtype=float), np.asarray(score1, dtype=float)]
    pos = y == 1
    neg = ~pos
    m, n = int(pos.sum()), int(neg.sum())
    if m < 2 or n < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan
    v10, v01, aucs = [], [], []
    for x in scores:
        xp, xn = x[pos], x[neg]
        phi = (xp[:, None] > xn[None, :]).astype(float) + 0.5 * (xp[:, None] == xn[None, :])
        v10.append(phi.mean(axis=1))
        v01.append(phi.mean(axis=0))
        aucs.append(float(phi.mean()))
    var = (np.cov(np.vstack(v10), ddof=1) / m + np.cov(np.vstack(v01), ddof=1) / n)
    variance = float(var[1, 1] + var[0, 0] - 2 * var[0, 1])
    delta = aucs[1] - aucs[0]
    if variance <= 0:
        return aucs[0], aucs[1], delta, np.nan, np.nan
    z = delta / np.sqrt(variance)
    return aucs[0], aucs[1], delta, float(z), float(2 * stats.norm.sf(abs(z)))


def incremental_validity(g):
    rows = []
    for outcome in SYMPTOMS:
        row = _ols_compare(g, outcome, GOLDS, GOLDS + ALPHAS)
        row.update({"analysis": "四金标准基础上加入四alpha", "effect": "delta_R2"})
        rows.append(row)
    for outcome in SYMPTOMS:
        sub = g[[outcome] + GOLDS + ALPHAS].dropna().copy()
        y = (sub[outcome] >= 10).astype(int)
        fits = []
        for cols in [GOLDS, GOLDS + ALPHAS]:
            fit = sm.GLM(y, _design(sub, cols), family=sm.families.Binomial()).fit()
            fits.append(fit)
        if not all(bool(fit.converged) for fit in fits):
            raise RuntimeError("临床结局Logistic模型未收敛：%s" % outcome)
        auc0, auc1, delta, z, p = _paired_delong(y, fits[0].predict(), fits[1].predict())
        r20 = 1 - fits[0].deviance / fits[0].null_deviance
        r21 = 1 - fits[1].deviance / fits[1].null_deviance
        rows.append({
            "outcome": outcome + "_clinical10", "n": int(len(sub)),
            "R2_base": float(r20), "R2_full": float(r21), "delta_R2": float(r21 - r20),
            "F": np.nan, "df_numerator": 4, "df_denominator": np.nan, "p": p,
            "AUC_base": auc0, "AUC_full": auc1, "delta_AUC": delta, "DeLong_z": z,
            "analysis": "四金标准基础上加入四alpha", "effect": "delta_AUC", "status": "ok",
        })
    return pd.DataFrame(rows)


def symptom_partial_correlations(g):
    rows = []
    for gold in GOLDS:
        for alpha in ALPHAS:
            sub = g[[gold, alpha] + SYMPTOMS].dropna()
            x0 = _design(sub, SYMPTOMS)
            rx = sm.OLS(sub[alpha].astype(float), x0).fit().resid
            ry = sm.OLS(sub[gold].astype(float), x0).fit().resid
            r = float(stats.pearsonr(rx, ry).statistic)
            # 两个变量均控制 SYMPTOMS 中的 k=2 个协变量；偏相关 t 检验须再扣除 k 个自由度。
            df = len(sub) - len(SYMPTOMS) - 2
            t = r * np.sqrt(df / (1 - r * r))
            p = float(2 * stats.t.sf(abs(t), df))
            rows.append({"criterion": gold, "alpha": alpha, "n": len(sub),
                         "partial_r": r, "t": t, "df": df, "p": p,
                         "controls": "PHQ-9+GAD-7"})
    out = pd.DataFrame(rows)
    out["p_bonferroni_16"] = np.minimum(1, out["p"] * len(out))
    return out


def _centered_profiles(d, cols):
    x = d[cols].astype(float)
    return x.sub(x.mean(axis=1), axis=0).add_prefix("profile_")


def _permutation_p(d, outcome, base_cols, profile_cols, observed_delta, seed=20261101, B=2000):
    """Freedman–Lane 残差置换：保留协变量与结局的关系，只检验新增剖面项。"""
    sub = d[[outcome] + base_cols + profile_cols].dropna()
    y = sub[outcome].to_numpy(dtype=float)
    x0 = _design(sub, base_cols).to_numpy(dtype=float)
    x1 = _design(sub, base_cols + profile_cols).to_numpy(dtype=float)
    fitted0 = x0 @ np.linalg.lstsq(x0, y, rcond=None)[0]
    residual0 = y - fitted0

    def delta(yy):
        r0 = yy - x0 @ np.linalg.lstsq(x0, yy, rcond=None)[0]
        r1 = yy - x1 @ np.linalg.lstsq(x1, yy, rcond=None)[0]
        return float(1 - (r1 @ r1) / (r0 @ r0))

    # delta() 是相对简化模型残差的增量，即 (R2_full-R2_base)/(1-R2_base)。
    # 将调用方传入的绝对 delta_R2 转到相同量纲，避免观测值与置换分布错配。
    observed_relative = delta(y)
    base_r2 = float(1 - (residual0 @ residual0) / (((y - y.mean()) @ (y - y.mean()))))
    expected_relative = observed_delta / (1 - base_r2)
    assert np.isclose(observed_relative, expected_relative, atol=1e-12, rtol=1e-10)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(B):
        permuted_y = fitted0 + rng.permutation(residual0)
        count += delta(permuted_y) >= observed_relative - 1e-12
    return (count + 1) / (B + 1), B


def ipsative_profiles(g):
    joint, coefficients, correlations = [], [], []
    for outcome in GOLDS:
        profile = _centered_profiles(g, ALPHAS)
        work = pd.concat([g[[outcome] + SYMPTOMS], profile], axis=1)
        full_profiles = list(profile.columns)
        # 四个居中剖面线性相关，省略其中一个即可张成同一三维空间。
        model_profiles = full_profiles[:-1]
        row = _ols_compare(work, outcome, SYMPTOMS, SYMPTOMS + model_profiles)
        perm_p, B = _permutation_p(work, outcome, SYMPTOMS, model_profiles, row["delta_R2"])
        row.update({"model": "full_4alpha_profile", "permutation_p": perm_p,
                    "permutations": B, "permutation_method": "Freedman-Lane",
                    "analysis": "组内四alpha居中剖面"})
        joint.append(row)

        reduced = _centered_profiles(g, ["alpha1", "alpha2", "alpha4"])
        work2 = pd.concat([g[[outcome] + SYMPTOMS], reduced], axis=1)
        reduced_profiles = list(reduced.columns)
        model2 = reduced_profiles[:-1]
        row2 = _ols_compare(work2, outcome, SYMPTOMS, SYMPTOMS + model2)
        perm2, _ = _permutation_p(work2, outcome, SYMPTOMS, model2, row2["delta_R2"])
        row2.update({"model": "exclude_alpha3_recentered", "permutation_p": perm2,
                     "permutations": B, "permutation_method": "Freedman-Lane",
                     "analysis": "排除alpha3后三alpha重新居中"})
        joint.append(row2)

        sub = work[[outcome] + SYMPTOMS + full_profiles].dropna()
        m0 = sm.OLS(sub[outcome], _design(sub, SYMPTOMS)).fit()
        m1 = sm.OLS(sub[outcome], _design(sub, SYMPTOMS + model_profiles)).fit()
        for term in model_profiles:
            coefficients.append({"criterion": outcome, "model": "full_4alpha_profile",
                                 "term": term, "beta": m1.params[term], "SE": m1.bse[term],
                                 "t": m1.tvalues[term], "p": m1.pvalues[term], "n": len(sub)})
        for term in full_profiles:
            r, p = stats.pearsonr(sub[term], sub[outcome])
            correlations.append({"criterion": outcome, "term": term, "n": len(sub), "r": r, "p": p})
    joint_df = pd.DataFrame(joint)
    joint_df["p_bonferroni_4criteria"] = np.minimum(1, joint_df["p"] * 4)
    corr_df = pd.DataFrame(correlations)
    corr_df["p_bonferroni_16"] = np.minimum(1, corr_df["p"] * len(corr_df))
    return joint_df, pd.DataFrame(coefficients), corr_df


def main():
    g = C.load_gold()
    assert g is not None and len(g) == 142
    C.save_tab(incremental_validity(g), "a41_Gold增量效度", "prolific", MODULE)
    C.save_tab(symptom_partial_correlations(g), "a41_症状偏相关_16对", "prolific", MODULE)
    joint, coef, corr = ipsative_profiles(g)
    C.save_tab(joint, "a41_组内剖面联合检验", "prolific", MODULE)
    C.save_tab(coef, "a41_组内剖面系数", "prolific", MODULE)
    C.save_tab(corr, "a41_组内剖面单相关", "prolific", MODULE)
    print("[a41] Gold增量效度、16对症状偏相关与组内剖面分析完成（n=142）")


if __name__ == "__main__":
    main()
