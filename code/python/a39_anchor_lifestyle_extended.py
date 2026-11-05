# a39_anchor_lifestyle_extended.py —— 锚定-算法偏差与生活方式扩展
# 全部为探索性分析；横断面路径不证明因果中介，不输出题目级专有编码。
import os

import numpy as np
import pandas as pd
from scipy.stats import f_oneway

import common as C

MODULE = "a39_anchor_lifestyle_extended"
ANCHOR = ["anchor1", "anchor2", "anchor3", "anchor4"]
ALPHA = ["alpha1", "alpha2", "alpha3", "alpha4"]
BOOT_B = int(os.environ.get("P3_BOOT_B", "2000"))
SEED = 3901


def _save(df, name, cohort):
    C.save_tab(df, name, cohort, MODULE)


def _status(cohort, analysis, message):
    return {"队列": cohort, "分析": analysis, "状态": message}


def _prepare(d):
    if any(c not in d for c in ANCHOR + ALPHA):
        return None
    out = d.copy()
    # 固定口径：anchor1–4 统一除以5至0–1；算法α沿用common官方基底。
    for col in ANCHOR:
        out[col + "_01"] = pd.to_numeric(out[col], errors="coerce") / 5.0
    anchored = [c + "_01" for c in ANCHOR]
    out["alpha_anch_mean"] = out[anchored].mean(axis=1, skipna=True)
    out["alpha_alg_mean"] = out[ALPHA].mean(axis=1, skipna=True)
    out["discrepancy"] = out["alpha_anch_mean"] - out["alpha_alg_mean"]
    out["abs_discrepancy"] = out["discrepancy"].abs()
    return out


def descriptive_reliability(d, cohort):
    anchored = [c + "_01" for c in ANCHOR]
    variables = anchored + ALPHA + ["alpha_anch_mean", "alpha_alg_mean", "discrepancy", "abs_discrepancy"]
    rows = []
    for col in variables:
        x = d[col].dropna()
        rows.append({"队列": cohort, "变量": col, "n": len(x), "均值": x.mean(),
                     "SD": x.std(ddof=1), "最小": x.min(), "P25": x.quantile(0.25),
                     "中位数": x.median(), "P75": x.quantile(0.75), "最大": x.max()})
    _save(pd.DataFrame(rows), "a39_锚定算法偏差_描述", cohort)
    reliability = []
    for label, cols in [("锚定四题_除5后", anchored), ("算法四alpha", ALPHA)]:
        complete = d[cols].dropna()
        reliability.append({"队列": cohort, "指标": label + " Cronbach_alpha",
                            "完整n": len(complete), "估计": C.cronbach_alpha(complete.values),
                            "状态": "已估计" if len(complete) >= 3 else "小样本"})
    _save(pd.DataFrame(reliability), "a39_Cronbach_alpha", cohort)


def correlations(d, cohort):
    rows = []
    pairs = [("alpha_anch_mean", "alpha_alg_mean", "锚定均值-算法均值"),
             ("discrepancy", "alpha_anch_mean", "偏差-锚定均值"),
             ("discrepancy", "alpha_alg_mean", "偏差-算法均值"),
             ("abs_discrepancy", "alpha_alg_mean", "绝对偏差-算法均值")]
    for a, b, label in pairs:
        sub = d[[a, b]].dropna()
        r, p = C.pearson_r(sub[a], sub[b]) if len(sub) >= 3 else (np.nan, np.nan)
        rows.append({"队列": cohort, "相关对": label, "变量1": a, "变量2": b,
                     "n": len(sub), "Pearson_r": r, "p": p,
                     "状态": "已估计" if len(sub) >= 3 else "小样本"})
    for a, b in zip(ANCHOR, ALPHA):
        ac = a + "_01"
        sub = d[[ac, b]].dropna()
        r, p = C.pearson_r(sub[ac], sub[b]) if len(sub) >= 3 else (np.nan, np.nan)
        rows.append({"队列": cohort, "相关对": "同维度锚定-算法", "变量1": ac, "变量2": b,
                     "n": len(sub), "Pearson_r": r, "p": p,
                     "状态": "已估计" if len(sub) >= 3 else "小样本"})
    _save(pd.DataFrame(rows), "a39_锚定算法_对应相关", cohort)


def discrepancy_outcomes(d, cohort):
    rows = []
    for outcome in ["phq9_total", "gad7_total", "ibs_positive", "ibs_symptom_count"]:
        if outcome not in d:
            continue
        for exposure in ["discrepancy", "abs_discrepancy"]:
            sub = d[[exposure, outcome]].dropna()
            r, p = C.pearson_r(sub[exposure], sub[outcome]) if len(sub) >= 3 else (np.nan, np.nan)
            rows.append({"队列": cohort, "偏差指标": exposure, "结局": outcome, "n": len(sub),
                         "Pearson_r": r, "p": p,
                         "说明": "二分类IBS时为点二列相关（Pearson等价）" if outcome == "ibs_positive" else "Pearson相关",
                         "状态": "已估计" if len(sub) >= 3 else "小样本"})
    _save(pd.DataFrame(rows), "a39_偏差与PHQ_GAD_IBS相关", cohort)


def discrepancy_anova(d, cohort):
    rows = []
    if "alpha_group" not in d:
        _save(pd.DataFrame([_status(cohort, "偏差alpha_group ANOVA", "alpha_group缺失")]),
              "a39_偏差_alpha_group_ANOVA", cohort)
        return
    for outcome in ["discrepancy", "abs_discrepancy", "alpha_anch_mean", "alpha_alg_mean"]:
        groups, arrays = [], []
        for group in ["耗竭型", "均衡型", "高意识型", "其他"]:
            values = d.loc[d["alpha_group"] == group, outcome].dropna()
            if len(values):
                groups.append(group)
                arrays.append(values)
        row = {"队列": cohort, "指标": outcome, "组数": len(groups),
               "各组n": ";".join("%s=%d" % (g, len(x)) for g, x in zip(groups, arrays)),
               "各组均值": ";".join("%s=%.6g" % (g, x.mean()) for g, x in zip(groups, arrays))}
        usable = [x for x in arrays if len(x) >= 2]
        if len(usable) >= 2:
            stat, p = f_oneway(*usable)
            row.update({"F": stat, "p": p, "n": sum(map(len, usable)), "状态": "已估计"})
        else:
            row["状态"] = "有效组不足或组内n<2，未检验"
        rows.append(row)
    _save(pd.DataFrame(rows), "a39_偏差_alpha_group_ANOVA", cohort)


def incremental_r2(d, cohort):
    import statsmodels.api as sm

    rows = []
    anchored = [c + "_01" for c in ANCHOR]
    for outcome in ["phq9_total", "gad7_total"]:
        if outcome not in d:
            continue
        sub = d[[outcome] + anchored + ALPHA].dropna()
        row = {"队列": cohort, "结局": outcome, "n": len(sub),
               "基线": "anchor1-4/5", "增量预测子": "alpha1-4"}
        if len(sub) < 20 or sub[outcome].std() == 0:
            row["状态"] = "小样本或结局零变异"
            rows.append(row)
            continue
        try:
            x0 = sm.add_constant(sub[anchored].astype(float), has_constant="add")
            x1 = sm.add_constant(sub[anchored + ALPHA].astype(float), has_constant="add")
            fit0 = sm.OLS(sub[outcome].astype(float), x0).fit()
            fit1 = sm.OLS(sub[outcome].astype(float), x1).fit()
            f_stat, p, df_diff = fit1.compare_f_test(fit0)
            row.update({"锚定基线R2": fit0.rsquared, "锚定基线调整R2": fit0.rsquared_adj,
                        "加入算法alpha后R2": fit1.rsquared, "完整模型调整R2": fit1.rsquared_adj,
                        "增量R2": fit1.rsquared - fit0.rsquared,
                        "增量F": f_stat, "增量df": df_diff, "增量p": p, "状态": "已估计"})
        except Exception as exc:
            row["状态"] = "模型失败：%s" % type(exc).__name__
        rows.append(row)
    _save(pd.DataFrame(rows), "a39_锚定基线到算法alpha_增量R2", cohort)


def sleep_quadratic(d, cohort):
    import statsmodels.api as sm

    rows = []
    covars = [c for c in ["smoke", "alcohol_binary"] if c in d]
    for outcome in ["phq9_total", "alpha3"]:
        if outcome not in d or "sleep_hour" not in d:
            continue
        sub = d[[outcome, "sleep_hour"] + covars].replace([np.inf, -np.inf], np.nan).dropna().copy()
        if len(sub) < 20 or sub["sleep_hour"].std() == 0 or sub[outcome].std() == 0:
            rows.append({"队列": cohort, "结局": outcome, "n": len(sub),
                         "状态": "小样本或睡眠/结局零变异"})
            continue
        center = sub["sleep_hour"].mean()
        sub["sleep_c"] = sub["sleep_hour"] - center
        sub["sleep_c2"] = sub["sleep_c"] ** 2
        for adjusted in [False, True]:
            if adjusted and not covars:
                continue
            preds = ["sleep_c", "sleep_c2"] + (covars if adjusted else [])
            X = sm.add_constant(sub[preds].astype(float), has_constant="add")
            try:
                fit = sm.OLS(sub[outcome].astype(float), X).fit()
                for pred in preds:
                    rows.append({"队列": cohort, "结局": outcome,
                                 "模型": "吸烟饮酒调整" if adjusted else "未调整",
                                 "预测子": pred, "coef": fit.params[pred], "SE": fit.bse[pred],
                                 "p": fit.pvalues[pred], "R2": fit.rsquared,
                                 "调整R2": fit.rsquared_adj, "n": len(sub), "睡眠中心化值": center,
                                 "二次项定义": "(sleep_hour-均值)^2", "状态": "已估计"})
            except Exception as exc:
                rows.append({"队列": cohort, "结局": outcome,
                             "模型": "吸烟饮酒调整" if adjusted else "未调整",
                             "n": len(sub), "状态": "模型失败：%s" % type(exc).__name__})
    _save(pd.DataFrame(rows), "a39_sleep_hour_二次项_PHQ_alpha3", cohort)


def _mediation_effects(values):
    # values = [alpha3, sleep_disorder_positive, PHQ]；二分类中介采用线性概率模型。
    X0 = np.column_stack([np.ones(len(values)), values[:, 0]])
    a = np.linalg.lstsq(X0, values[:, 1], rcond=None)[0][1]
    total = np.linalg.lstsq(X0, values[:, 2], rcond=None)[0][1]
    X1 = np.column_stack([np.ones(len(values)), values[:, 0], values[:, 1]])
    direct_b = np.linalg.lstsq(X1, values[:, 2], rcond=None)[0]
    direct, b = direct_b[1], direct_b[2]
    return np.array([a, b, a * b, total, direct])


def sleep_mediation(d, cohort):
    rows = []
    if cohort == "jianshu":
        rows.append(_status(cohort, "alpha3→sleep→PHQ中介", "不适用：见数睡眠阳性按任务要求作废"))
        _save(pd.DataFrame(rows), "a39_Prolific_sleep中介_bootstrap", cohort)
        return
    needed = ["alpha3", "sleep_disorder_positive", "phq9_total"]
    if any(c not in d for c in needed):
        rows.append(_status(cohort, "alpha3→sleep→PHQ中介", "Prolific所需变量缺失，不适用"))
        _save(pd.DataFrame(rows), "a39_Prolific_sleep中介_bootstrap", cohort)
        return
    sub = d[needed].dropna()
    if len(sub) < 30 or sub["sleep_disorder_positive"].nunique() < 2 or sub["alpha3"].std() == 0:
        rows.append(_status(cohort, "alpha3→sleep→PHQ中介", "n<30、睡眠阳性无变异或alpha3零变异，不适用"))
        _save(pd.DataFrame(rows), "a39_Prolific_sleep中介_bootstrap", cohort)
        return
    values = sub[needed].values.astype(float)
    estimates = _mediation_effects(values)
    rng = np.random.default_rng(SEED)
    boots = []
    for _ in range(BOOT_B):
        idx = rng.integers(0, len(values), len(values))
        try:
            boots.append(_mediation_effects(values[idx]))
        except np.linalg.LinAlgError:
            pass
    boots = np.asarray(boots)
    labels = ["a_alpha3到sleep_线性概率", "b_sleep到PHQ_控制alpha3",
              "间接效应a乘b", "总效应", "直接效应"]
    for j, label in enumerate(labels):
        if len(boots) >= 100:
            lo, hi = np.quantile(boots[:, j], [0.025, 0.975])
            p = min(1.0, 2 * min(np.mean(boots[:, j] <= 0), np.mean(boots[:, j] >= 0)))
        else:
            lo, hi, p = np.nan, np.nan, np.nan
        rows.append({"队列": cohort, "效应": label, "估计": estimates[j],
                     "bootstrap_CI低": lo, "bootstrap_CI高": hi,
                     "bootstrap_p_双侧符号": p, "有效bootstrap_n": len(boots), "n": len(sub),
                     "方法边界": "二分类睡眠中介用线性概率模型，PHQ用OLS；探索性横断面乘积法，不作因果确认"})
    _save(pd.DataFrame(rows), "a39_Prolific_sleep中介_bootstrap", cohort)


def discrepancy_covariate_sensitivity(d, cohort):
    import statsmodels.api as sm

    rows = []
    covars = [c for c in ["smoke", "alcohol_binary"] if c in d]
    if not covars:
        _save(pd.DataFrame([_status(cohort, "吸烟饮酒协变量敏感性", "smoke/alcohol_binary均缺失")]),
              "a39_偏差结局_吸烟饮酒敏感性", cohort)
        return
    for outcome in ["phq9_total", "gad7_total", "ibs_positive"]:
        if outcome not in d:
            continue
        sub = d[["discrepancy", outcome] + covars].replace([np.inf, -np.inf], np.nan).dropna()
        binary = outcome == "ibs_positive"
        if len(sub) < 20 or sub[outcome].std() == 0 or (binary and sub[outcome].nunique() < 2):
            rows.append({"队列": cohort, "结局": outcome, "n": len(sub),
                         "状态": "小样本或结局无变异"})
            continue
        for adjusted in [False, True]:
            preds = ["discrepancy"] + (covars if adjusted else [])
            X = sm.add_constant(sub[preds].astype(float), has_constant="add")
            try:
                if binary:
                    fit = sm.GLM(sub[outcome].astype(float), X, family=sm.families.Binomial()).fit()
                    r2 = 1 - fit.llf / fit.llnull
                else:
                    fit = sm.OLS(sub[outcome].astype(float), X).fit()
                    r2 = fit.rsquared
                rows.append({"队列": cohort, "结局": outcome,
                             "模型": "吸烟饮酒调整" if adjusted else "未调整",
                             "偏差coef": fit.params["discrepancy"], "SE": fit.bse["discrepancy"],
                             "OR_仅二分类": np.exp(fit.params["discrepancy"]) if binary else np.nan,
                             "p": fit.pvalues["discrepancy"],
                             "R2或McFadden_R2": r2, "n": len(sub),
                             "协变量": "+".join(covars) if adjusted else "无", "状态": "已估计；两模型同一样本"})
            except Exception as exc:
                rows.append({"队列": cohort, "结局": outcome,
                             "模型": "吸烟饮酒调整" if adjusted else "未调整", "n": len(sub),
                             "状态": "模型失败：%s" % type(exc).__name__})
    _save(pd.DataFrame(rows), "a39_偏差结局_吸烟饮酒敏感性", cohort)


def smoking_alcohol_joint_groups(d, cohort):
    if "smoke" not in d or "alcohol_binary" not in d:
        _save(pd.DataFrame([_status(cohort, "吸烟饮酒联合组", "smoke或alcohol_binary缺失")]),
              "a39_吸烟饮酒联合组", cohort)
        return
    work = d.copy()
    work["联合组"] = np.where(work["smoke"].isna() | work["alcohol_binary"].isna(), np.nan,
                         np.where((work["smoke"] > 0) & (work["alcohol_binary"] == 1), "吸烟+饮酒",
                         np.where((work["smoke"] > 0) & (work["alcohol_binary"] == 0), "吸烟+不饮酒",
                         np.where((work["smoke"] <= 0) & (work["alcohol_binary"] == 1), "不吸烟+饮酒", "不吸烟+不饮酒"))))
    rows, tests = [], []
    metrics = ["discrepancy", "abs_discrepancy", "phq9_total", "gad7_total", "alpha3", "sleep_hour"]
    for group in ["不吸烟+不饮酒", "不吸烟+饮酒", "吸烟+不饮酒", "吸烟+饮酒"]:
        sub = work[work["联合组"] == group]
        row = {"队列": cohort, "联合组": group, "n_分组定义完整": len(sub)}
        for metric in metrics:
            if metric in sub:
                x = sub[metric].dropna()
                row[metric + "_n"] = len(x)
                row[metric + "_均值"] = x.mean() if len(x) else np.nan
                row[metric + "_SD"] = x.std(ddof=1) if len(x) > 1 else np.nan
        rows.append(row)
    for metric in metrics:
        if metric not in work:
            continue
        arrays = [work.loc[work["联合组"] == g, metric].dropna()
                  for g in ["不吸烟+不饮酒", "不吸烟+饮酒", "吸烟+不饮酒", "吸烟+饮酒"]]
        usable = [x for x in arrays if len(x) >= 2]
        row = {"队列": cohort, "指标": metric, "有效组数": len(usable),
               "n": sum(map(len, usable))}
        if len(usable) >= 2:
            stat, p = f_oneway(*usable)
            row.update({"ANOVA_F": stat, "p": p, "状态": "已估计"})
        else:
            row["状态"] = "有效组不足，未检验"
        tests.append(row)
    _save(pd.DataFrame(rows), "a39_吸烟饮酒联合组", cohort)
    _save(pd.DataFrame(tests), "a39_吸烟饮酒联合组_ANOVA", cohort)


def main():
    notes = pd.DataFrame([
        ["推断属性", "全部结果为探索性；除非另有可核查预注册，不得表述为确认性证据。"],
        ["锚定统一", "anchor1–4固定除以5至0–1；算法α沿用common官方基底；discrepancy=锚定均值-算法均值。"],
        ["中介限制", "仅Prolific在sleep_disorder_positive可用且有变异时估计；见数睡眠阳性作废并标不适用。二分类中介用线性概率乘积法。"],
        ["因果边界", "横断面路径、增量R²与bootstrap区间不证明时间顺序、因果中介或临床效用。"],
        ["联合组", "吸烟=smoke>0；饮酒=alcohol_binary==1，构成四个互斥组；缺失者不强行归类。"],
        ["隐私与资产", "仅输出聚合统计，不输出个人记录、题目级专有编码或映射。"]
    ], columns=["项目", "说明"])
    for cohort in C.COHORTS:
        raw = C.load_cohort(cohort)
        if raw is None:
            _save(pd.DataFrame([_status(cohort, "数据加载", "基底未就绪；未运行且不伪造结果")]),
                  "a39_运行状态", cohort)
            continue
        d = _prepare(raw)
        if d is None:
            _save(pd.DataFrame([_status(cohort, "变量准备", "锚定或alpha变量缺失，未伪造结果")]),
                  "a39_运行状态", cohort)
            continue
        statuses = []
        jobs = [("描述与Cronbach", lambda: descriptive_reliability(d, cohort)),
                ("锚定算法相关", lambda: correlations(d, cohort)),
                ("偏差与结局相关", lambda: discrepancy_outcomes(d, cohort)),
                ("alpha_group ANOVA", lambda: discrepancy_anova(d, cohort)),
                ("增量R2", lambda: incremental_r2(d, cohort)),
                ("睡眠二次项", lambda: sleep_quadratic(d, cohort)),
                ("Prolific睡眠中介", lambda: sleep_mediation(d, cohort)),
                ("吸烟饮酒敏感性", lambda: discrepancy_covariate_sensitivity(d, cohort)),
                ("吸烟饮酒联合组", lambda: smoking_alcohol_joint_groups(d, cohort))]
        for name, job in jobs:
            try:
                job()
                statuses.append(_status(cohort, name, "完成；具体不可估计项见各表状态"))
            except Exception as exc:
                statuses.append(_status(cohort, name, "跳过：%s：%s" % (type(exc).__name__, str(exc)[:160])))
        _save(pd.DataFrame(statuses), "a39_运行状态", cohort)
        _save(notes.assign(队列=cohort), "a39_分析口径与边界", cohort)
        print("[a39] %s 锚定与生活方式扩展完成（探索性）" % cohort)


if __name__ == "__main__":
    main()
