# a25_replication_extension.py —— 复现核定与扩展
# 保留公开 H_trend = α1² 的输出级内部一致性、肥胖 OR、BMI 相关、
# 数值锚点信度与下游 profile 描述统计；不读取内部字段或重构后台。
import numpy as np
import pandas as pd
import common as C

MODULE = "a25_replication_extension"


def htrend_internal_consistency(d, cohort):
    """公开 H_trend = α1² 与后台 H 的横断面关联及替代模型 AIC。"""
    import statsmodels.api as sm
    sub = d[["H_score", "alpha1", "H_trend"]].dropna()
    a1sq = sub["alpha1"].to_numpy(dtype=float) ** 2
    h = sub["H_score"].to_numpy(dtype=float)
    mq = sm.OLS(h, sm.add_constant(a1sq, has_constant="add")).fit()
    ml = sm.OLS(h, sm.add_constant(sub["alpha1"].to_numpy(dtype=float), has_constant="add")).fit()
    r = float(np.corrcoef(a1sq, h)[0, 1])
    return pd.DataFrame([{
        "cohort": cohort, "n": len(sub),
        "r_alpha1sq_vs_H": round(r, 4), "R2_alpha1sq_vs_H": round(r ** 2, 4),
        "AIC_quadratic": round(float(mq.aic), 2),
        "AIC_linear": round(float(ml.aic), 2),
        "note": "H_trend is the public downstream transformation alpha1^2; identity itself is not evidence",
    }])


def bmi_incremental(d, cohort):
    """年龄、性别、教育基线上逐加入一个α；输出H12/SM3.3完整聚合表。"""
    import statsmodels.api as sm
    rows = []
    covars = ["age_group", "gender", "education"]
    for alpha in ["alpha1", "alpha2", "alpha3", "alpha4"]:
        sub = d[["bmi", alpha] + covars].dropna().copy()
        x0 = pd.get_dummies(sub[covars], columns=covars, drop_first=True, dtype=float)
        x0 = sm.add_constant(x0, has_constant="add")
        x1 = x0.copy()
        x1[alpha] = sub[alpha].astype(float)
        m0 = sm.OLS(sub["bmi"].astype(float), x0).fit()
        m1 = sm.OLS(sub["bmi"].astype(float), x1).fit()
        rows.append({"队列": cohort, "alpha": alpha, "n": len(sub),
                     "beta": m1.params[alpha], "SE": m1.bse[alpha],
                     "t": m1.tvalues[alpha], "p": m1.pvalues[alpha],
                     "delta_R2": m1.rsquared-m0.rsquared, "base_R2": m0.rsquared})
    return pd.DataFrame(rows)


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        rows = []
        # 肥胖 OR（H_score 预测 obese）
        if "obese" in d.columns:
            sub = d[["obese", "H_score"]].dropna()
            if sub["obese"].nunique() >= 2 and len(sub) > 30:
                ors = C.logit_or(sub, "obese", ["H_score"])
                rows.append(("肥胖 OR(每单位H)", round(ors["H_score"][0], 4),
                             round(ors["H_score"][1], 5)))
        # BMI 相关
        if "bmi" in d.columns:
            r, pv = C.pearson_r(d["bmi"], d["H_score"])
            rows.append(("BMI 与 H 相关 r", round(r, 4), round(pv, 6)))
        # 锚定信度：仅在受控 L3 含四个锚点题时计算；公开 L2 不伪造空结果。
        anchor_cols = ["anchor1", "anchor2", "anchor3", "anchor4"]
        if all(c in d.columns for c in anchor_cols):
            anc = d[anchor_cols].dropna()
            if len(anc) > 3:
                rows.append(("锚定四题 Cronbach α", round(C.cronbach_alpha(anc.values), 4), np.nan))
        # top-1 profile
        top = d["profile"].value_counts().head(3)
        for k, v in top.items():
            rows.append(("最高频profile=%s" % k, int(v), np.nan))
        tab = pd.DataFrame(rows, columns=["指标", "数值", "p"])
        tab.insert(0, "队列", cohort)
        C.save_tab(tab, "a25_复现扩展", cohort, MODULE)
        C.save_tab(htrend_internal_consistency(d, cohort), "T25_htrend_replication", cohort, MODULE)
        C.save_tab(bmi_incremental(d, cohort), "a25_BMI四alpha增量回归", cohort, MODULE)
        print("[a25] %s 完成" % cohort)
    print("[a25] 复现核定与扩展完成")


if __name__ == "__main__":
    main()
