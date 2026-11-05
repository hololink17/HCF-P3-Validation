# a04_predictive.py —— 预测效度 / 判别力 AUC（移植 R 10_analysis_04_predictive.R + 34_disease_or.R + 36_seven_factors.R）
# 每个因子（α1-4 / carrier a-d / H_score）对多个二分类结局的判别力 AUC（DeLong CI）。
# 同时输出七疾病 OR（逻辑回归，每病单独模型）。两队列。
import numpy as np
import pandas as pd
import common as C

MODULE = "a04_predictive"
PREDICTORS = ["alpha1", "alpha2", "alpha3", "alpha4",
              "carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score"]
OUTCOMES = ["depression", "ibs_positive", "chronic_pain", "heart_disease",
            "hypertension", "diabetes", "hyperlipidemia", "obese"]


def auc_table(d, cohort):
    rows = []
    for outcome in OUTCOMES:
        if outcome not in d.columns:
            continue
        y = d[outcome].dropna()
        if y.nunique() < 2:
            continue
        for pred in PREDICTORS:
            s = d[pred].dropna()
            yy = y.reindex(s.index) if False else d.loc[s.index, outcome]
            m = yy.notna() & s.notna()
            if m.sum() < 20:
                continue
            auc, lo, hi = C.discrimination(yy[m].values, s[m].values)
            rows.append((outcome, pred, round(auc, 4), round(lo, 4), round(hi, 4),
                        int(m.sum()), int(yy[m].sum())))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "AUC", "CI低", "CI高", "n", "阳性数"])
    out.insert(0, "队列", cohort)
    return out


def disease_or(d, cohort):
    rows = []
    for outcome in OUTCOMES:
        if outcome not in d.columns:
            continue
        y = d[outcome].dropna()
        if y.nunique() < 2:
            continue
        for pred in PREDICTORS:
            sub = d[[outcome, pred]].dropna()
            if sub[outcome].nunique() < 2 or len(sub) < 30:
                continue
            ors = C.logit_or(sub, outcome, [pred])
            if pred in ors:
                o, p = ors[pred]
                rows.append((outcome, pred, round(o, 4), round(p, 5),
                             int((sub[outcome] == 1).sum()), len(sub)))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "OR", "p", "阳性数", "n"])
    out.insert(0, "队列", cohort)
    return out


def univariate_per_sd(d, cohort):
    """九读数×全部结局：健康方向AUC及单变量per-SD Logistic OR/CI。"""
    import statsmodels.api as sm
    rows = []
    for outcome in OUTCOMES:
        for pred in PREDICTORS:
            sub = d[[outcome, pred]].dropna().copy()
            if len(sub) < 30 or sub[outcome].nunique() < 2 or sub[pred].std(ddof=1) == 0:
                continue
            y = sub[outcome].astype(float)
            z = (sub[pred]-sub[pred].mean())/sub[pred].std(ddof=1)
            x = sm.add_constant(z.astype(float), has_constant="add")
            try:
                fit = sm.GLM(y, x, family=sm.families.Binomial()).fit()
                if not bool(fit.converged):
                    status = "nonconverged"
                    rows.append({"队列": cohort, "结局": outcome, "预测因子": pred,
                                 "n": len(sub), "状态": status})
                    continue
                auc, lo, hi = C.discrimination(y, sub[pred])
                ci = fit.conf_int().loc[pred]
                rows.append({"队列": cohort, "结局": outcome, "预测因子": pred,
                             "AUC": auc, "AUC_CI低": lo, "AUC_CI高": hi,
                             "OR_per_SD": np.exp(fit.params[pred]),
                             "OR_CI低": np.exp(ci[0]), "OR_CI高": np.exp(ci[1]),
                             "p": fit.pvalues[pred], "n": len(sub),
                             "阳性数": int(y.sum()), "状态": "ok"})
            except Exception as e:
                rows.append({"队列": cohort, "结局": outcome, "预测因子": pred,
                             "n": len(sub), "状态": type(e).__name__})
    return pd.DataFrame(rows)


def adjusted_disease_or(d, cohort):
    """年龄、性别、教育校正的OR及95%CI；与R34同口径，失败单列不中断。"""
    import statsmodels.api as sm
    rows = []
    covars = [c for c in ["age_group", "gender", "education"] if c in d.columns]
    for outcome in OUTCOMES:
        if outcome not in d.columns:
            continue
        for pred in PREDICTORS:
            cols = [outcome, pred] + covars
            sub = d[cols].dropna().copy()
            if len(sub) < 40 or sub[outcome].nunique() < 2:
                continue
            try:
                X = pd.get_dummies(sub[[pred] + covars], columns=covars,
                                   drop_first=True, dtype=float)
                X = sm.add_constant(X.astype(float), has_constant="add")
                # 同一二项最大似然模型使用 GLM/IRLS；避免 Newton 数值容差造成假未收敛。
                fit = sm.GLM(sub[outcome].astype(float), X,
                             family=sm.families.Binomial()).fit()
                if not bool(fit.converged):
                    rows.append((outcome, pred, np.nan, np.nan, np.nan, np.nan,
                                 np.nan, np.nan, np.nan, len(sub),
                                 int(sub[outcome].sum()), "nonconverged"))
                    continue
                ci = fit.conf_int().loc[pred]
                rows.append((outcome, pred, float(np.exp(fit.params[pred])),
                             float(np.exp(ci[0])), float(np.exp(ci[1])),
                             float(np.exp(0.1 * fit.params[pred])),
                             float(np.exp(0.1 * ci[0])), float(np.exp(0.1 * ci[1])),
                             float(fit.pvalues[pred]), len(sub), int(sub[outcome].sum()), "ok"))
            except Exception as e:
                rows.append((outcome, pred, np.nan, np.nan, np.nan, np.nan,
                             np.nan, np.nan, np.nan, len(sub),
                             int(sub[outcome].sum()), type(e).__name__))
    out = pd.DataFrame(rows, columns=["结局", "预测因子", "校正OR", "CI低", "CI高",
                                      "校正OR_per_0_1", "CI低_per_0_1", "CI高_per_0_1",
                                      "p", "n", "阳性数", "状态"])
    # 明确两种检验族，避免正文引用口径不明：四α×8结局，以及九读数×8结局。
    from statsmodels.stats.multitest import multipletests
    estimable = out["p"].notna()
    out["q_BH_9readouts_8outcomes"] = np.nan
    out.loc[estimable, "q_BH_9readouts_8outcomes"] = multipletests(
        out.loc[estimable, "p"], method="fdr_bh")[1]
    alpha_rows = estimable & out["预测因子"].isin(["alpha1", "alpha2", "alpha3", "alpha4"])
    out["q_BH_4alpha_8outcomes"] = np.nan
    out.loc[alpha_rows, "q_BH_4alpha_8outcomes"] = multipletests(
        out.loc[alpha_rows, "p"], method="fdr_bh")[1]
    out.insert(0, "队列", cohort)
    return out


def symptom_adjusted_disease(d, cohort):
    """年龄+性别+教育+PHQ-9+GAD-7基线上，逐一加入九读数的增量模型。"""
    import statsmodels.api as sm
    from scipy import stats
    rows = []
    covars = ["age_group", "gender", "education"]
    symptoms = ["phq9_total", "gad7_total"]
    for outcome in ["ibs_positive", "chronic_pain", "heart_disease",
                    "hypertension", "diabetes", "hyperlipidemia", "obese"]:
        for pred in PREDICTORS:
            cols = [outcome, pred] + covars + symptoms
            sub = d[cols].dropna().copy()
            y = sub[outcome].astype(float)
            x = pd.get_dummies(sub[covars], columns=covars, drop_first=True, dtype=float)
            for symptom in symptoms:
                x[symptom] = sub[symptom].astype(float)
            x0 = sm.add_constant(x.astype(float), has_constant="add")
            x1 = x0.copy()
            x1[pred] = sub[pred].astype(float)
            try:
                m0 = sm.GLM(y, x0, family=sm.families.Binomial()).fit()
                m1 = sm.GLM(y, x1, family=sm.families.Binomial()).fit()
                if not (bool(m0.converged) and bool(m1.converged)):
                    status = "nonconverged"
                    rows.append({"队列": cohort, "结局": outcome, "预测因子": pred,
                                 "n": len(sub), "状态": status})
                    continue
                lr = max(0.0, 2*(m1.llf-m0.llf))
                lr_p = float(stats.chi2.sf(lr, 1))
                auc0 = C.auc_delong(y, m0.predict())[0]
                auc1 = C.auc_delong(y, m1.predict())[0]
                ci = m1.conf_int().loc[pred]
                rows.append({"队列": cohort, "结局": outcome, "预测因子": pred, "n": len(sub),
                             "AUC_base": auc0, "AUC_full": auc1, "delta_AUC": auc1-auc0,
                             "LR_chi2": lr, "LR_p": lr_p, "coef_p": m1.pvalues[pred],
                             "OR_per_0_1": np.exp(.1*m1.params[pred]),
                             "CI低_per_0_1": np.exp(.1*ci[0]),
                             "CI高_per_0_1": np.exp(.1*ci[1]), "状态": "ok"})
            except Exception as e:
                rows.append({"队列": cohort, "结局": outcome, "预测因子": pred,
                             "n": len(sub), "状态": type(e).__name__})
    out = pd.DataFrame(rows)
    from statsmodels.stats.multitest import multipletests
    out["q_BH_within_outcome_9"] = np.nan
    for _, idx in out.groupby("结局").groups.items():
        idx = list(idx)
        estimable = out.loc[idx, "LR_p"].notna()
        if estimable.any():
            chosen = out.loc[idx].loc[estimable].index
            out.loc[chosen, "q_BH_within_outcome_9"] = multipletests(
                out.loc[chosen, "LR_p"], method="fdr_bh")[1]
    return out


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        C.save_tab(auc_table(d, cohort), "a04_判别力AUC", cohort, MODULE)
        C.save_tab(disease_or(d, cohort), "a04_疾病OR", cohort, MODULE)
        C.save_tab(univariate_per_sd(d, cohort), "a04_九读数疾病perSD", cohort, MODULE)
        C.save_tab(adjusted_disease_or(d, cohort), "a04_疾病OR_校正", cohort, MODULE)
        C.save_tab(symptom_adjusted_disease(d, cohort), "a04_控症状后疾病增量", cohort, MODULE)
    print("[a04] 预测效度 AUC + 疾病 OR + 控症状增量完成")


if __name__ == "__main__":
    main()
