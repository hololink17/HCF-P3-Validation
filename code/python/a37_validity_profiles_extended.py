# a37_validity_profiles_extended.py —— α 结构、理论型谱与增量效度扩展
# 全部为探索性分析；仅使用 common.py 已生成的聚合分析变量，不重算专有编码。
import itertools
import os

import numpy as np
import pandas as pd
from scipy.stats import chi2, chi2_contingency, f_oneway

import common as C

MODULE = "a37_validity_profiles_extended"
ALPHA = ["alpha1", "alpha2", "alpha3", "alpha4"]
BOOT_B = int(os.environ.get("P3_BOOT_B", "2000"))
SEED = 3701


def _status(cohort, analysis, message):
    return {"队列": cohort, "分析": analysis, "状态": message}


def _save(df, name, cohort):
    C.save_tab(df, name, cohort, MODULE)


def _add_fdr(store, cohort, family, test, p, n=np.nan):
    if p is not None and np.isfinite(p):
        store.append({"队列": cohort, "检验族": family, "检验": test,
                      "p": float(p), "n": n})


def _bh(store):
    if not store:
        return pd.DataFrame(columns=["队列", "检验族", "检验", "p", "n", "BH_q", "BH拒绝"])
    out = pd.DataFrame(store)
    out["BH_q"] = np.nan
    out["BH拒绝"] = False
    # 在每一队列内对所有主要探索检验统一校正，而非按分析块分别校正。
    for cohort, idx in out.groupby("队列").groups.items():
        pos = out.index.get_indexer(idx)
        p = out.loc[idx, "p"].astype(float).values
        order = np.argsort(p)
        ranked = p[order]
        m = len(ranked)
        q = ranked * m / np.arange(1, m + 1)
        q = np.minimum.accumulate(q[::-1])[::-1]
        q_full = np.empty(m)
        q_full[order] = np.minimum(q, 1.0)
        out.loc[idx, "BH_q"] = q_full
        out.loc[idx, "BH拒绝"] = q_full < 0.05
    return out.sort_values(["队列", "BH_q", "p"])


def simplex(d, cohort, fdr):
    pair_rows = []
    for i, j in itertools.combinations(range(4), 2):
        a, b = ALPHA[i], ALPHA[j]
        if a not in d or b not in d:
            continue
        sub = d[[a, b]].dropna()
        if len(sub) < 10:
            pair_rows.append({"队列": cohort, "因子1": a, "因子2": b, "间隔": j - i,
                              "类型": "相邻" if j - i == 1 else "非相邻", "n": len(sub),
                              "r": np.nan, "p": np.nan, "期望方向": "r>0",
                              "方向符合": np.nan, "状态": "小样本/缺失，未检验"})
            continue
        r, p = C.pearson_r(sub[a], sub[b])
        pair_rows.append({"队列": cohort, "因子1": a, "因子2": b, "间隔": j - i,
                          "类型": "相邻" if j - i == 1 else "非相邻", "n": len(sub),
                          "r": r, "p": p, "期望方向": "r>0；总体上相邻相关强于非相邻",
                          "方向符合": bool(r > 0), "状态": "已估计"})
        _add_fdr(fdr, cohort, "simplex相关", "%s-%s" % (a, b), p, len(sub))
    pairs = pd.DataFrame(pair_rows)
    _save(pairs, "a37_simplex_相关对", cohort)
    if pairs.empty or pairs["r"].dropna().empty:
        summary = pd.DataFrame([_status(cohort, "simplex", "无可用相关对")])
    else:
        adj = pairs.loc[pairs["类型"] == "相邻", "r"].dropna()
        non = pairs.loc[pairs["类型"] == "非相邻", "r"].dropna()
        diff = adj.mean() - non.mean() if len(adj) and len(non) else np.nan
        summary = pd.DataFrame([{
            "队列": cohort, "相邻对数": len(adj), "非相邻对数": len(non),
            "相邻平均r": adj.mean() if len(adj) else np.nan,
            "非相邻平均r": non.mean() if len(non) else np.nan,
            "差值_相邻减非相邻": diff,
            "明确期望": "全部相关为正，且相邻平均r-非相邻平均r>0",
            "总体判据符合": bool(diff > 0) if np.isfinite(diff) else np.nan,
            "说明": "该均值差为描述性 simplex 判据，不等同于单调递减的逐对证明。"
        }])
    _save(summary, "a37_simplex_摘要", cohort)


def symptom_regressions(d, cohort, fdr):
    import statsmodels.api as sm

    rows = []
    for outcome in ["phq9_total", "gad7_total"]:
        if outcome not in d:
            continue
        specs = [("单独", [a]) for a in ALPHA] + [("四alpha联合", ALPHA)]
        for model_type, preds in specs:
            if any(p not in d for p in preds):
                continue
            sub = d[[outcome] + preds].dropna()
            if len(sub) < max(12, len(preds) + 5) or sub[outcome].std() == 0:
                rows.append({"队列": cohort, "结局": outcome, "模型": model_type,
                             "预测子": "+".join(preds), "n": len(sub), "状态": "小样本/零变异"})
                continue
            X = sm.add_constant(sub[preds].astype(float), has_constant="add")
            res = sm.OLS(sub[outcome].astype(float), X).fit()
            for pred in preds:
                rows.append({"队列": cohort, "结局": outcome, "模型": model_type,
                             "预测子": pred, "n": len(sub), "coef": res.params[pred],
                             "SE": res.bse[pred], "p": res.pvalues[pred],
                             "R2": res.rsquared, "调整R2": res.rsquared_adj,
                             "模型F_p": res.f_pvalue, "期望方向": "负向（α越高，症状越低）",
                             "状态": "已估计"})
                _add_fdr(fdr, cohort, "PHQ/GAD回归", "%s~%s:%s" %
                         (outcome, model_type, pred), res.pvalues[pred], len(sub))
    _save(pd.DataFrame(rows), "a37_PHQ_GAD_alpha回归", cohort)


def profile_counts(d, cohort):
    if "profile" not in d:
        _save(pd.DataFrame([_status(cohort, "81型谱", "profile缺失")]),
              "a37_81理论profile计数", cohort)
        return
    levels = ["低", "中", "高"]
    theoretical = ["-".join(x) for x in itertools.product(levels, repeat=4)]
    counts = d["profile"].value_counts()
    out = pd.DataFrame({"profile": theoretical})
    out["n"] = out["profile"].map(counts).fillna(0).astype(int)
    out["占完整profile样本pct"] = 100 * out["n"] / max(int(out["n"].sum()), 1)
    out["实际出现"] = out["n"] > 0
    out.insert(0, "队列", cohort)
    out["理论profile总数"] = 81
    out["完整profile_n"] = int(out["n"].sum())
    out["profile缺失_n"] = int(d["profile"].isna().sum())
    _save(out, "a37_81理论profile计数", cohort)


def _demographics(d):
    parts = []
    for col in ["age_group", "gender"]:
        if col in d:
            z = pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=float)
            parts.append(z)
    return pd.concat(parts, axis=1) if parts else pd.DataFrame(index=d.index)


def interaction_lr(d, cohort, fdr):
    import statsmodels.api as sm

    rows = []
    for outcome in ["ibs_positive", "chronic_pain"]:
        if outcome not in d or any(a not in d for a in ["alpha2", "alpha3"]):
            continue
        demo = _demographics(d)
        work = pd.concat([d[[outcome, "alpha2", "alpha3"]], demo], axis=1)
        work = work.replace([np.inf, -np.inf], np.nan).dropna()
        if len(work) < 30 or work[outcome].nunique() < 2:
            rows.append({"队列": cohort, "结局": outcome, "n": len(work),
                         "状态": "缺失/小样本/结局无病例-对照变异，LR不适用"})
            continue
        a2 = work["alpha2"] - work["alpha2"].mean()
        a3 = work["alpha3"] - work["alpha3"].mean()
        base = pd.DataFrame({"alpha2_c": a2, "alpha3_c": a3}, index=work.index)
        demo_cols = [c for c in work.columns if c.startswith("age_group_") or c.startswith("gender_")]
        base = pd.concat([base, work[demo_cols].astype(float)], axis=1)
        full = base.copy()
        full["alpha2_x_alpha3"] = a2 * a3
        X0 = sm.add_constant(base, has_constant="add")
        X1 = sm.add_constant(full, has_constant="add")
        try:
            fit0 = sm.GLM(work[outcome].astype(float), X0, family=sm.families.Binomial()).fit()
            fit1 = sm.GLM(work[outcome].astype(float), X1, family=sm.families.Binomial()).fit()
            lr = max(float(fit0.deviance - fit1.deviance), 0.0)
            df_lr = int(X1.shape[1] - X0.shape[1])
            p_lr = chi2.sf(lr, df_lr)
            term = "alpha2_x_alpha3"
            rows.append({"队列": cohort, "结局": outcome, "n": len(work),
                         "交互coef": fit1.params[term], "交互SE": fit1.bse[term],
                         "交互OR": np.exp(fit1.params[term]), "Wald_p": fit1.pvalues[term],
                         "LR_chi2": lr, "LR_df": df_lr, "LR_p": p_lr,
                         "全模型McFadden_R2": 1 - fit1.llf / fit1.llnull,
                         "协变量": "age_group+gender（虚拟变量）；α2/α3已均值中心化",
                         "状态": "已估计"})
            _add_fdr(fdr, cohort, "α2×α3交互LR", outcome, p_lr, len(work))
        except Exception as exc:
            rows.append({"队列": cohort, "结局": outcome, "n": len(work),
                         "状态": "模型未收敛/不可估计：%s" % type(exc).__name__})
    _save(pd.DataFrame(rows), "a37_alpha2xalpha3_交互LR", cohort)


def ibs_subset_auc(d, cohort):
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    rows = []
    if "ibs_positive" not in d or any(a not in d for a in ALPHA):
        _save(pd.DataFrame([_status(cohort, "15组合AUC", "IBS或α变量缺失")]),
              "a37_IBS_15个alpha子集_AUC_CV", cohort)
        return
    full = d[["ibs_positive"] + ALPHA].dropna()
    y = full["ibs_positive"].astype(int)
    for k in range(1, 5):
        for subset in itertools.combinations(ALPHA, k):
            label = "+".join(subset)
            row = {"队列": cohort, "alpha子集": label, "预测子数": k, "n": len(full),
                   "阳性n": int(y.sum()) if len(y) else 0, "设定CV折数": 5}
            if len(full) < 20 or y.nunique() < 2:
                row["状态"] = "小样本或IBS无病例-对照变异，AUC不可恢复"
                rows.append(row)
                continue
            X = full[list(subset)].astype(float)
            model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))
            try:
                model.fit(X, y)
                risk = model.predict_proba(X)[:, 1]
                row["表观AUC"] = roc_auc_score(y, risk)
                if len(full) >= 30 and y.value_counts().min() >= 5:
                    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
                    cv_risk = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
                    row["5折CV_AUC"] = roc_auc_score(y, cv_risk)
                    row["状态"] = "已估计"
                else:
                    row["状态"] = "表观AUC已估计；每类不足5例或n<30，5折CV不适用"
            except Exception as exc:
                row["状态"] = "模型失败：%s" % type(exc).__name__
            rows.append(row)
    _save(pd.DataFrame(rows), "a37_IBS_15个alpha子集_AUC_CV", cohort)


def _chain_coefficients(x):
    # x 列为 α1、α2、α3、α4；链式特定间接效应用于 a12*a23*a34。
    def slope(y, predictors):
        X = np.column_stack([np.ones(len(x))] + [x[:, j] for j in predictors])
        return np.linalg.lstsq(X, y, rcond=None)[0]

    a = slope(x[:, 1], [0])[1]
    b = slope(x[:, 2], [0, 1])[2]
    c = slope(x[:, 3], [0, 1, 2])[3]
    return a * b * c


def chain_path(d, cohort, fdr):
    import statsmodels.api as sm

    rows = []
    if any(a not in d for a in ALPHA):
        _save(pd.DataFrame([_status(cohort, "链式路径", "α变量缺失")]),
              "a37_alpha链式路径_bootstrap", cohort)
        return
    sub = d[ALPHA].dropna()
    if len(sub) < 30 or (sub.std() == 0).any():
        _save(pd.DataFrame([_status(cohort, "链式路径", "n<30或α零变异")]),
              "a37_alpha链式路径_bootstrap", cohort)
        return
    z = (sub - sub.mean()) / sub.std(ddof=0)
    equations = [("alpha2", ["alpha1"]), ("alpha3", ["alpha1", "alpha2"]),
                 ("alpha4", ["alpha1", "alpha2", "alpha3"])]
    fits = {}
    for outcome, preds in equations:
        fit = sm.OLS(z[outcome], sm.add_constant(z[preds], has_constant="add")).fit()
        fits[outcome] = fit
        for pred in preds:
            rows.append({"队列": cohort, "效应": "%s→%s" % (pred, outcome),
                         "类型": "标准化路径（控制前序变量）", "估计": fit.params[pred],
                         "SE": fit.bse[pred], "p": fit.pvalues[pred], "n": len(z)})
            _add_fdr(fdr, cohort, "链式路径", "%s→%s" % (pred, outcome),
                     fit.pvalues[pred], len(z))
    a = fits["alpha2"].params["alpha1"]
    b = fits["alpha3"].params["alpha2"]
    c = fits["alpha4"].params["alpha3"]
    indirect = a * b * c
    rng = np.random.default_rng(SEED)
    values = z[ALPHA].values
    boots = []
    for _ in range(BOOT_B):
        idx = rng.integers(0, len(values), len(values))
        try:
            boots.append(_chain_coefficients(values[idx]))
        except np.linalg.LinAlgError:
            pass
    boots = np.asarray(boots, dtype=float)
    boots = boots[np.isfinite(boots)]
    lo, hi = (np.quantile(boots, [0.025, 0.975]) if len(boots) >= 100 else (np.nan, np.nan))
    p_boot = (2 * min(np.mean(boots <= 0), np.mean(boots >= 0))
              if len(boots) else np.nan)
    rows.append({"队列": cohort, "效应": "α1→α2→α3→α4", "类型": "链式特定间接效应a*b*c",
                 "估计": indirect, "bootstrap_CI低": lo, "bootstrap_CI高": hi,
                 "bootstrap_p_双侧符号": min(p_boot, 1.0) if np.isfinite(p_boot) else np.nan,
                 "有效bootstrap_n": len(boots), "n": len(z),
                 "状态": "Python statsmodels/numpy线性路径与百分位bootstrap"})
    _add_fdr(fdr, cohort, "链式路径", "α1→α2→α3→α4间接效应", p_boot, len(z))
    _save(pd.DataFrame(rows), "a37_alpha链式路径_bootstrap", cohort)


def _clinical_grade(value, scale):
    if pd.isna(value):
        return np.nan
    if value < 5:
        return "0-4_最轻"
    if value < 10:
        return "5-9_轻度"
    if value < 15:
        return "10-14_中度"
    if scale == "PHQ" and value < 20:
        return "15-19_中重度"
    return ("20-27_重度" if scale == "PHQ" else "15-21_重度")


def clinical_crosstabs(d, cohort, fdr):
    rows = []
    for scale, col in [("PHQ", "phq9_total"), ("GAD", "gad7_total")]:
        if col not in d or "alpha_group" not in d:
            continue
        work = d[[col, "alpha_group"]].dropna().copy()
        work["临床分级"] = work[col].map(lambda x: _clinical_grade(x, scale))
        if work.empty:
            continue
        order = (["0-4_最轻", "5-9_轻度", "10-14_中度", "15-19_中重度", "20-27_重度"]
                 if scale == "PHQ" else ["0-4_最轻", "5-9_轻度", "10-14_中度", "15-21_重度"])
        tab = pd.crosstab(work["alpha_group"], work["临床分级"]).reindex(columns=order, fill_value=0)
        for group in tab.index:
            for grade in order:
                rows.append({"队列": cohort, "量表": scale, "alpha组合": group,
                             "临床分级": grade, "n": int(tab.loc[group, grade]),
                             "组内pct": 100 * tab.loc[group, grade] / max(tab.loc[group].sum(), 1)})
        if tab.shape[0] >= 2 and tab.shape[1] >= 2 and tab.values.sum() > 0:
            try:
                stat, p, dof, _ = chi2_contingency(tab)
                rows.append({"队列": cohort, "量表": scale, "alpha组合": "整体卡方",
                             "临床分级": "α组合×临床分级", "n": int(tab.values.sum()),
                             "chi2": stat, "df": dof, "p": p})
                _add_fdr(fdr, cohort, "α组合×临床分级", scale, p, int(tab.values.sum()))
            except ValueError:
                pass
    _save(pd.DataFrame(rows), "a37_alpha组合xPHQ_GAD临床分级", cohort)


def three_group_tests(d, cohort, fdr):
    rows = []
    groups = ["潜在抑郁", "典型抑郁", "真正健康"]
    if "group_type" not in d:
        _save(pd.DataFrame([_status(cohort, "三组检验", "group_type缺失")]),
              "a37_潜在典型真正健康_三组检验", cohort)
        return
    work = d[d["group_type"].isin(groups)].copy()
    for col in ["phq9_total", "gad7_total", "ibs_symptom_count", "alpha1", "alpha2", "alpha3", "alpha4"]:
        if col not in work:
            continue
        arrays = [work.loc[work["group_type"] == g, col].dropna() for g in groups]
        row = {"队列": cohort, "指标": col, "检验": "一元ANOVA",
               **{"n_%s" % g: len(a) for g, a in zip(groups, arrays)},
               **{"均值_%s" % g: a.mean() if len(a) else np.nan for g, a in zip(groups, arrays)}}
        if min(map(len, arrays)) >= 2:
            stat, p = f_oneway(*arrays)
            row.update({"统计量": stat, "p": p, "n": sum(map(len, arrays)), "状态": "已估计"})
            _add_fdr(fdr, cohort, "潜在/典型/真正健康三组", col, p, row["n"])
        else:
            row["状态"] = "至少一组n<2，未检验"
        rows.append(row)
    for col in ["ibs_positive", "chronic_pain"]:
        if col not in work:
            continue
        tab = pd.crosstab(work["group_type"], work[col]).reindex(index=groups, columns=[0, 1], fill_value=0)
        row = {"队列": cohort, "指标": col, "检验": "Pearson卡方", "n": int(tab.values.sum()),
               **{"n_%s" % g: int(tab.loc[g].sum()) for g in groups},
               **{"阳性率_%s" % g: (tab.loc[g, 1] / tab.loc[g].sum()
                                  if tab.loc[g].sum() else np.nan) for g in groups}}
        if tab.values.sum() and (tab.sum(axis=1) > 0).sum() >= 2:
            stat, p, dof, _ = chi2_contingency(tab)
            row.update({"统计量": stat, "df": dof, "p": p, "状态": "已估计；稀疏单元格时仅作探索"})
            _add_fdr(fdr, cohort, "潜在/典型/真正健康三组", col, p, row["n"])
        else:
            row["状态"] = "有效组不足，未检验"
        rows.append(row)
    _save(pd.DataFrame(rows), "a37_潜在典型真正健康_三组检验", cohort)


def main():
    notes = pd.DataFrame([
        ["推断属性", "全部结果为探索性；除非另有可核查预注册，不得表述为确认性证据。"],
        ["simplex期望", "α相关预期为正，且相邻相关的总体强度高于非相邻相关。"],
        ["链式实现边界", "仅以Python statsmodels/numpy估计线性路径及百分位bootstrap间接效应；不是lavaan SEM，未提供也不声称CFI/RMSEA。"],
        ["链式解释边界", "链式方向来自理论设定，横断面关联与bootstrap区间不证明时间顺序或因果中介。"],
        ["CV", "IBS组合模型使用固定5折分层CV；类别数不足时不恢复或伪造CV。"],
        ["多重检验", "主要探索检验在每一队列内统一BH-FDR；描述性AUC、型谱计数与净重抽样不另造p值。"],
        ["隐私与资产", "仅输出聚合统计，不输出个人记录、题目级专有编码或映射。"]
    ], columns=["项目", "说明"])
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        if d is None:
            _save(pd.DataFrame([_status(cohort, "数据加载", "基底未就绪；未运行且不伪造结果")]),
                  "a37_运行状态", cohort)
            continue
        fdr, statuses = [], []
        jobs = [
            ("simplex", lambda: simplex(d, cohort, fdr)),
            ("PHQ/GAD回归", lambda: symptom_regressions(d, cohort, fdr)),
            ("81型谱", lambda: profile_counts(d, cohort)),
            ("α2×α3交互", lambda: interaction_lr(d, cohort, fdr)),
            ("15组合AUC/CV", lambda: ibs_subset_auc(d, cohort)),
            ("链式路径", lambda: chain_path(d, cohort, fdr)),
            ("临床分级", lambda: clinical_crosstabs(d, cohort, fdr)),
            ("三组检验", lambda: three_group_tests(d, cohort, fdr)),
        ]
        for name, job in jobs:
            try:
                job()
                statuses.append(_status(cohort, name, "完成；具体不可估计项见各表状态"))
            except Exception as exc:
                statuses.append(_status(cohort, name, "跳过：%s：%s" % (type(exc).__name__, str(exc)[:160])))
        _save(_bh(fdr), "a37_主要探索检验_BH_FDR", cohort)
        _save(pd.DataFrame(statuses), "a37_运行状态", cohort)
        _save(notes.assign(队列=cohort), "a37_分析口径与边界", cohort)
        print("[a37] %s 扩展分析完成（探索性）" % cohort)


if __name__ == "__main__":
    main()
