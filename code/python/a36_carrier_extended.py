# a36_carrier_extended.py —— 四载体扩展分析（全部探索性，非确认性结论）
import warnings

import numpy as np
import pandas as pd
from scipy.stats import f_oneway
import statsmodels.api as sm
from statsmodels.stats.multicomp import pairwise_tukeyhsd

import common as C

MODULE = "a36_carrier_extended"
CARRIERS = ["carrier_a", "carrier_b", "carrier_c", "carrier_d"]
DISEASES = [
    "ibs_positive", "chronic_pain", "heart_disease",
    "hypertension", "diabetes", "hyperlipidemia",
]
DEMOGRAPHICS = ["age_group", "gender", "education"]
LEVELS = ["低", "中", "高"]
MIN_N = 40
MIN_OUTCOME_CELL = 5
MIN_EVENTS_PER_PARAMETER = 5
EXPLORATORY = "exploratory"


def _with_meta(df, cohort):
    out = df.copy()
    out.insert(0, "cohort", cohort)
    out["exploratory"] = EXPLORATORY
    return out


def _save(df, name, cohort):
    C.save_tab(_with_meta(df, cohort), name, cohort, MODULE)


def _error(errors, analysis, reason, outcome="", predictor="", model=""):
    errors.append({
        "analysis": analysis, "outcome": outcome, "predictor": predictor,
        "model": model, "reason": reason,
    })


def _prepare(d, errors):
    d = d.copy()
    numeric = set(CARRIERS + DISEASES + [
        "sleep_hour", "bmi", "obese", "smoke",
        "alcohol_freq", "alcohol_binary",
    ])
    for col in numeric:
        if col in d.columns:
            d[col] = pd.to_numeric(d[col], errors="coerce")
    cutoffs = []
    for carrier in CARRIERS:
        cat_col = carrier + "_cat"
        if carrier not in d.columns or d[carrier].dropna().nunique() < 2:
            d[cat_col] = np.nan
            _error(errors, "载体P33/P67分档", "载体缺失或无变异", predictor=carrier)
            cutoffs.append({"carrier": carrier, "P33": np.nan, "P67": np.nan,
                            "n_nonmissing": 0})
            continue
        s = d[carrier]
        p33, p67 = np.nanquantile(s.values, [1 / 3, 2 / 3])
        cat = pd.Series(np.nan, index=d.index, dtype=object)
        cat[s.notna() & (s <= p33)] = "低"
        cat[s.notna() & (s > p33) & (s < p67)] = "中"
        cat[s.notna() & (s >= p67)] = "高"
        d[cat_col] = cat
        if p33 == p67:
            _error(errors, "载体P33/P67分档", "P33=P67，部分档位可能为空",
                   predictor=carrier, model="分档保留但后续小样本检查")
        cutoffs.append({"carrier": carrier, "P33": p33, "P67": p67,
                        "n_nonmissing": int(s.notna().sum())})
    return d, pd.DataFrame(cutoffs)


def _reduced_design(X, required):
    """依次去除零方差与完全共线项；核心载体若不可估计则明确报错。"""
    keep = []
    dropped = 0
    const = np.ones((len(X), 1))
    rank = 1
    for col in X.columns:
        values = X[col].astype(float).values
        if np.nanstd(values) == 0:
            if col in required:
                raise ValueError("核心载体零方差，无法估计")
            dropped += 1
            continue
        trial = keep + [col]
        matrix = np.column_stack([const, X[trial].astype(float).values])
        new_rank = np.linalg.matrix_rank(matrix)
        if new_rank == rank + 1:
            keep.append(col)
            rank = new_rank
        else:
            if col in required:
                raise ValueError("核心载体与截距或其他预测子完全共线")
            dropped += 1
    return X[keep], dropped


def _fit_logistic(d, outcome, carriers, adjusted, model_label):
    """成功时返回每个核心载体的OR；失败时抛出异常，由调用方写入错误表。"""
    covars = [c for c in DEMOGRAPHICS if adjusted and c in d.columns]
    needed = [outcome] + carriers + covars
    missing = [c for c in needed if c not in d.columns]
    if missing:
        raise ValueError("缺少模型列：" + ",".join(missing))
    dd = d[needed].dropna().copy()
    dd = dd[dd[outcome].isin([0, 1])]
    if len(dd) < MIN_N:
        raise ValueError("完整病例n<%d" % MIN_N)
    y = dd[outcome].astype(float)
    events = int((y == 1).sum())
    non_events = int((y == 0).sum())
    if min(events, non_events) < MIN_OUTCOME_CELL:
        raise ValueError("结局较小组<%d，避免不稳定Logistic估计" % MIN_OUTCOME_CELL)

    X = dd[carriers].astype(float).copy()
    if covars:
        dummy = pd.get_dummies(dd[covars].astype(str), drop_first=True, dtype=float)
        X = pd.concat([X, dummy], axis=1)
    if len(dd) <= X.shape[1] + 2:
        raise ValueError("相对参数数目样本量不足")
    X, dropped = _reduced_design(X, set(carriers))
    X = sm.add_constant(X, has_constant="add")
    if min(events, non_events) / X.shape[1] < MIN_EVENTS_PER_PARAMETER:
        raise ValueError("较小结局组/参数数<%d，模型样本量不足" % MIN_EVENTS_PER_PARAMETER)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = sm.Logit(y, X).fit(disp=0, maxiter=200)
    if not result.mle_retvals.get("converged", True):
        raise ValueError("Logistic未收敛")
    params = result.params
    bse = result.bse
    if not np.isfinite(params).all() or not np.isfinite(bse).all():
        raise ValueError("模型系数或标准误非有限值")
    if (np.abs(params.drop("const", errors="ignore")) > 25).any() or \
            (bse.drop("const", errors="ignore") > 10).any():
        raise ValueError("可能存在完全/准完全分离或稀疏类别，估计不稳定")
    if any(w.category.__name__ in {"PerfectSeparationWarning", "PerfectSeparationError"}
           for w in caught):
        raise ValueError("检测到完全分离")

    ci = result.conf_int()
    rows = []
    for carrier in carriers:
        if carrier not in params:
            raise ValueError("核心载体未能进入最终模型：" + carrier)
        rows.append({
            "outcome": outcome,
            "predictor": carrier,
            "model": model_label,
            "OR_per_1_unit": float(np.exp(params[carrier])),
            "CI_low": float(np.exp(ci.loc[carrier, 0])),
            "CI_high": float(np.exp(ci.loc[carrier, 1])),
            "p": float(result.pvalues[carrier]),
            "n": int(len(dd)),
            "events": events,
            "non_events": non_events,
            "n_parameters": int(X.shape[1]),
            "events_per_parameter": events / X.shape[1],
            "dropped_collinear_or_constant_covariates": dropped,
            "inference": "exploratory_not_confirmatory",
        })
    return rows


def _logistic_tables(d, errors):
    rows = []
    for outcome in DISEASES:
        if outcome not in d.columns:
            _error(errors, "单载体Logistic", "结局列缺失", outcome=outcome)
            continue
        for carrier in CARRIERS:
            for adjusted, label in [(False, "crude"), (True, "adjusted")]:
                try:
                    rows.extend(_fit_logistic(d, outcome, [carrier], adjusted, label))
                except Exception as exc:
                    _error(errors, "单载体Logistic", type(exc).__name__ + ": " + str(exc),
                           outcome=outcome, predictor=carrier, model=label)
        for adjusted, label in [(False, "joint_crude"), (True, "joint_adjusted")]:
            try:
                rows.extend(_fit_logistic(d, outcome, CARRIERS, adjusted, label))
            except Exception as exc:
                _error(errors, "四载体联合Logistic", type(exc).__name__ + ": " + str(exc),
                       outcome=outcome, predictor="carrier_a-d", model=label)
    return pd.DataFrame(rows)


def _profile_table(d):
    cats = [c + "_cat" for c in CARRIERS]
    complete = d[cats].notna().all(axis=1)
    profile = d.loc[complete, cats].agg("-".join, axis=1)
    out = profile.value_counts().rename_axis("carrier_profile").reset_index(name="n")
    if len(out):
        out["pct_of_complete_profiles"] = 100 * out["n"] / out["n"].sum()
    out["n_complete_profiles"] = int(complete.sum())
    out["n_total"] = len(d)
    return out


def _correlation_rows(d, carrier, outcomes, function_label, errors):
    rows = []
    for outcome in outcomes:
        if carrier not in d.columns or outcome not in d.columns:
            _error(errors, function_label + "相关", "载体或结局列缺失",
                   outcome=outcome, predictor=carrier)
            continue
        pair = d[[carrier, outcome]].dropna()
        if len(pair) < 3 or pair[carrier].nunique() < 2 or pair[outcome].nunique() < 2:
            _error(errors, function_label + "相关", "完整病例不足或变量无变异",
                   outcome=outcome, predictor=carrier)
            continue
        for method, fun in [("Pearson", C.pearson_r), ("Spearman", C.spearman_r)]:
            try:
                estimate, p = fun(pair[carrier].values, pair[outcome].values)
                rows.append({
                    "function": function_label, "carrier": carrier, "outcome": outcome,
                    "method": method, "n": len(pair), "estimate": estimate, "p": p,
                })
            except Exception as exc:
                _error(errors, function_label + "相关_" + method,
                       type(exc).__name__ + ": " + str(exc),
                       outcome=outcome, predictor=carrier)
    return rows


def _group_analyses(d, carrier, outcomes, function_label, errors):
    cat_col = carrier + "_cat"
    descriptive_rows, anova_rows, tukey_rows = [], [], []
    for outcome in outcomes:
        if outcome not in d.columns:
            _error(errors, function_label + "组间比较", "结局列缺失",
                   outcome=outcome, predictor=carrier)
            continue
        sub = d[[cat_col, outcome]].dropna()
        groups = []
        for level in LEVELS:
            values = sub.loc[sub[cat_col] == level, outcome].astype(float)
            groups.append(values)
            descriptive_rows.append({
                "function": function_label, "carrier": carrier, "outcome": outcome,
                "carrier_level": level, "n": len(values),
                "mean_or_rate": values.mean() if len(values) else np.nan,
                "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
                "small_sample_n_lt_5": bool(len(values) < 5),
            })
        try:
            if any(len(g) < 2 for g in groups):
                raise ValueError("低/中/高至少一组n<2")
            if sub[outcome].nunique() < 2:
                raise ValueError("结局无变异")
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                test = f_oneway(*groups)
            if not np.isfinite(test.statistic) or not np.isfinite(test.pvalue):
                raise ValueError("ANOVA估计非有限值")
            anova_rows.append({
                "function": function_label, "carrier": carrier, "outcome": outcome,
                "F": test.statistic, "p": test.pvalue, "n": int(len(sub)),
                "note": "obese为0/1率的探索性ANOVA" if outcome == "obese" else "",
            })
            tukey = pairwise_tukeyhsd(sub[outcome].astype(float).values,
                                      sub[cat_col].values, alpha=0.05)
            for values in tukey.summary().data[1:]:
                tukey_rows.append({
                    "function": function_label, "carrier": carrier, "outcome": outcome,
                    "group_1": values[0], "group_2": values[1],
                    "mean_difference": values[2], "p_adjusted": values[3],
                    "CI_low": values[4], "CI_high": values[5], "reject": values[6],
                })
        except Exception as exc:
            _error(errors, function_label + "ANOVA/Tukey",
                   type(exc).__name__ + ": " + str(exc), outcome=outcome, predictor=carrier)
    return descriptive_rows, anova_rows, tukey_rows


def _stratified_table(d, alcohol_col):
    outcomes = [c for c in ["sleep_hour", "smoke", alcohol_col] if c and c in d.columns]
    rows = []
    if "alpha_group" not in d.columns:
        return pd.DataFrame(rows)
    for carrier in CARRIERS:
        cat_col = carrier + "_cat"
        for level in LEVELS:
            for alpha_group in sorted(d["alpha_group"].dropna().astype(str).unique()):
                sub = d[(d[cat_col] == level) & (d["alpha_group"].astype(str) == alpha_group)]
                for outcome in outcomes:
                    values = sub[outcome].dropna().astype(float)
                    rows.append({
                        "carrier": carrier, "carrier_level": level,
                        "alpha_group": alpha_group, "outcome": outcome,
                        "n": len(values),
                        "mean_or_rate": values.mean() if len(values) else np.nan,
                        "sd": values.std(ddof=1) if len(values) > 1 else np.nan,
                        "small_sample_n_lt_5": bool(len(values) < 5),
                        "model": "carrier_x_alpha_group_stratified_descriptive",
                    })
    return pd.DataFrame(rows)


def _run_cohort(d, cohort):
    errors = []
    d, cutoffs = _prepare(d, errors)
    alcohol_col = "alcohol_freq" if "alcohol_freq" in d.columns else (
        "alcohol_binary" if "alcohol_binary" in d.columns else "")

    logistic = _logistic_tables(d, errors)
    profile = _profile_table(d)
    stratified = _stratified_table(d, alcohol_col)

    _save(cutoffs, "a36_各队列载体P33P67切点", cohort)
    _save(logistic, "a36_疾病Logistic_OR", cohort)
    _save(profile, "a36_载体组合profile分布", cohort)
    _save(stratified, "a36_载体乘alpha_group分层表", cohort)
    note = pd.DataFrame([
        {"item": "推断属性", "note": "所有结果均为exploratory；显著性不作为确认性结论或因果证据"},
        {"item": "载体分档", "note": "每个队列分别按本队列P33/P67划为低/中/高"},
        {"item": "Logistic", "note": "OR为连续载体每增加1单位；校正模型将age_group、gender、education作为分类哑变量"},
        {"item": "联合模型", "note": "joint_crude与joint_adjusted同时纳入carrier_a-d；完全共线、未收敛、分离或小样本模型写入错误表"},
        {"item": "公开边界", "note": "不指定或恢复任何carrier与内部语义、题目、症状或功能之间的映射；仅进行中性读数的下游统计"},
        {"item": "分层探索", "note": "载体×alpha_group采用分层聚合表，不作交互显著性的确认性解释；alcohol使用" + (alcohol_col or "不可用")},
    ])
    _save(note, "a36_分析说明", cohort)
    _save(pd.DataFrame(errors, columns=["analysis", "outcome", "predictor", "model", "reason"]),
          "a36_分析错误", cohort)


def main():
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        if d is None:
            _save(pd.DataFrame([{"analysis": "A36", "outcome": "", "predictor": "",
                                 "model": "", "reason": "分析基底未就绪；未伪造运行结果"}]),
                  "a36_分析错误", cohort)
            continue
        _run_cohort(d, cohort)
        print("[a36] %s 四载体扩展探索分析完成" % cohort)
    print("[a36] 双队列处理完成；所有推断均标记为exploratory")


if __name__ == "__main__":
    main()
