# a38_ibs_clinical_extended.py —— IBS 人群归因、判别与临床效用扩展
# 全部为探索性分析；PAR/DCA/NNT 不构成治疗获益或临床推荐。
import os

import numpy as np
import pandas as pd

import common as C

MODULE = "a38_ibs_clinical_extended"
BOOT_B = int(os.environ.get("P3_BOOT_B", "2000"))
SEED = 3801
THRESHOLDS = np.round(np.arange(0.01, 0.501, 0.01), 2)
MAJOR_COMORBIDITY = ["chronic_pain", "heart_disease", "hypertension", "diabetes", "hyperlipidemia"]


def _save(df, name, cohort):
    C.save_tab(df, name, cohort, MODULE)


def _status(cohort, analysis, message):
    return {"队列": cohort, "分析": analysis, "状态": message}


def _demographics(d):
    parts = []
    for col in ["age_group", "gender"]:
        if col in d:
            parts.append(pd.get_dummies(d[col], prefix=col, drop_first=True, dtype=float))
    return pd.concat(parts, axis=1) if parts else pd.DataFrame(index=d.index)


def _prepare(d, adjusted):
    if "ibs_positive" not in d or "alpha3" not in d:
        return None, None, "IBS或alpha3缺失"
    demo = _demographics(d) if adjusted else pd.DataFrame(index=d.index)
    work = pd.concat([d[["ibs_positive", "alpha3"]], demo], axis=1)
    work = work.replace([np.inf, -np.inf], np.nan).dropna()
    if len(work) < 20 or work["ibs_positive"].nunique() < 2:
        return None, None, "小样本或IBS无病例-对照变异，ROC不可恢复"
    predictors = ["alpha3"] + [c for c in work.columns if c.startswith("age_group_") or c.startswith("gender_")]
    X = work[predictors].astype(float)
    y = work["ibs_positive"].astype(int)
    return X, y, "已估计"


def _new_model():
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000))


def _cv_risk(model, X, y):
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    if len(y) < 30 or y.value_counts().min() < 5:
        return None, "每类不足5例或n<30，固定5折CV不适用"
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED)
    pred = cross_val_predict(model, X, y, cv=cv, method="predict_proba")[:, 1]
    return pred, "已估计"


def _bootstrap_auc(y, risk, seed):
    from sklearn.metrics import roc_auc_score
    rng = np.random.default_rng(seed)
    y = np.asarray(y, dtype=int)
    risk = np.asarray(risk, dtype=float)
    vals = []
    for _ in range(BOOT_B):
        idx = rng.integers(0, len(y), len(y))
        if len(np.unique(y[idx])) == 2:
            vals.append(roc_auc_score(y[idx], risk[idx]))
    if len(vals) < 100:
        return np.nan, np.nan, len(vals)
    lo, hi = np.quantile(vals, [0.025, 0.975])
    return lo, hi, len(vals)


def roc_models(d, cohort):
    from sklearn.metrics import roc_auc_score, roc_curve

    summaries, curves = [], []
    for adjusted, label in [(False, "alpha3单独"), (True, "alpha3+年龄性别")]:
        X, y, state = _prepare(d, adjusted)
        base = {"队列": cohort, "模型": label, "设定CV折数": 5}
        if X is None:
            summaries.append({**base, "状态": state})
            continue
        try:
            model = _new_model()
            model.fit(X, y)
            risk = model.predict_proba(X)[:, 1]
            cv_risk, cv_state = _cv_risk(_new_model(), X, y)
            lo, hi, nb = _bootstrap_auc(y, risk, SEED + int(adjusted))
            summaries.append({**base, "n": len(y), "阳性n": int(y.sum()),
                              "表观AUC": roc_auc_score(y, risk),
                              "AUC_bootstrap_CI低": lo, "AUC_bootstrap_CI高": hi,
                              "有效bootstrap_n": nb,
                              "5折CV_AUC": roc_auc_score(y, cv_risk) if cv_risk is not None else np.nan,
                              "CV状态": cv_state, "状态": "已估计；bootstrap为表观评分百分位区间"})
            fpr, tpr, thresholds = roc_curve(y, risk)
            for i in range(len(fpr)):
                curves.append({"队列": cohort, "模型": label, "n": len(y),
                               "FPR": fpr[i], "TPR": tpr[i],
                               "风险阈值": thresholds[i] if np.isfinite(thresholds[i]) else np.nan})
        except Exception as exc:
            summaries.append({**base, "n": len(y), "状态": "模型失败：%s" % type(exc).__name__})
    _save(pd.DataFrame(summaries), "a38_alpha3_ROC_CV_bootstrap", cohort)
    _save(pd.DataFrame(curves), "a38_alpha3_ROC坐标", cohort)


def population_attributable_fraction(d, cohort):
    row = {"队列": cohort, "暴露固定定义": "暴露=alpha3>=本队列P67；非暴露=alpha3<P67",
           "PAR公式": "(P总-P非暴露)/P总 = 1-P非暴露/P总；高α3为保护性暴露时PAR可为负"}
    if "ibs_positive" not in d or "alpha3" not in d:
        row["状态"] = "IBS或alpha3缺失"
        _save(pd.DataFrame([row]), "a38_IBS_PAR_alpha3P67", cohort)
        return
    sub = d[["ibs_positive", "alpha3"]].dropna().copy()
    if len(sub) < 20:
        row.update({"n": len(sub), "状态": "小样本，PAR不适用"})
        _save(pd.DataFrame([row]), "a38_IBS_PAR_alpha3P67", cohort)
        return
    cut = sub["alpha3"].quantile(2 / 3)
    exposed = (sub["alpha3"] >= cut).astype(int).values
    y = sub["ibs_positive"].astype(int).values

    def calc(yy, ee):
        if len(yy) == 0 or (ee == 0).sum() == 0 or yy.mean() == 0:
            return np.nan, np.nan, np.nan, np.nan, np.nan
        p_total = yy.mean()
        p_exp = yy[ee == 1].mean() if (ee == 1).any() else np.nan
        p_unexp = yy[ee == 0].mean()
        par = (p_total - p_unexp) / p_total
        rr = p_exp / p_unexp if p_unexp > 0 else np.nan
        return p_total, p_exp, p_unexp, par if np.isfinite(par) else np.nan, rr

    # 返回值依次为总风险、暴露风险、非暴露风险、PAR、RR。
    result = calc(y, exposed)
    p_total, p_exp, p_unexp, par, rr = result
    rng = np.random.default_rng(SEED + 10)
    boots = []
    for _ in range(BOOT_B):
        idx = rng.integers(0, len(y), len(y))
        value = calc(y[idx], exposed[idx])[3]
        if np.isfinite(value):
            boots.append(value)
    lo, hi = (np.quantile(boots, [0.025, 0.975]) if len(boots) >= 100 else (np.nan, np.nan))
    row.update({"n": len(sub), "alpha3_P67": cut, "暴露n": int(exposed.sum()),
                "非暴露n": int((1 - exposed).sum()), "IBS总风险": p_total,
                "暴露组IBS风险": p_exp, "非暴露组IBS风险": p_unexp,
                "RR_暴露vs非暴露": rr, "PAR": par,
                "PAR_bootstrap_CI低": lo, "PAR_bootstrap_CI高": hi,
                "有效bootstrap_n": len(boots), "状态": "已估计" if np.isfinite(par) else "PAR不可估计"})
    _save(pd.DataFrame([row]), "a38_IBS_PAR_alpha3P67", cohort)


def sensitivity(d, cohort):
    from sklearn.metrics import roc_auc_score

    rows = []
    masks = {"全体": pd.Series(True, index=d.index)}
    if "med_psych" in d:
        masks["排除精神用药"] = d["med_psych"].fillna(0) != 1
    comorb = pd.Series(False, index=d.index)
    available = [c for c in MAJOR_COMORBIDITY if c in d]
    for col in available:
        comorb = comorb | (d[col] == 1).fillna(False)
    if available:
        masks["排除主要共病"] = ~comorb
    if "med_psych" in d and available:
        masks["同时排除精神用药与主要共病"] = (d["med_psych"].fillna(0) != 1) & ~comorb
    for scenario, mask in masks.items():
        ds = d.loc[mask]
        for adjusted, label in [(False, "alpha3单独"), (True, "alpha3+年龄性别")]:
            X, y, state = _prepare(ds, adjusted)
            row = {"队列": cohort, "敏感性样本": scenario, "模型": label,
                   "排除定义": "精神用药=med_psych==1；主要共病=" + "+".join(available)}
            if X is None:
                row["状态"] = state  # 无病例-对照变异时不伪造亚组ROC。
                rows.append(row)
                continue
            try:
                model = _new_model().fit(X, y)
                risk = model.predict_proba(X)[:, 1]
                auc, lo, hi = C.auc_delong(y, risk)
                row.update({"n": len(y), "阳性n": int(y.sum()), "表观AUC": auc,
                            "DeLong_CI低": lo, "DeLong_CI高": hi, "状态": "已估计"})
            except Exception as exc:
                row["状态"] = "模型失败：%s" % type(exc).__name__
            rows.append(row)
    _save(pd.DataFrame(rows), "a38_IBS_排除用药共病敏感性", cohort)


def _net_benefit(y, risk, thresholds):
    y = np.asarray(y, dtype=int)
    risk = np.asarray(risk, dtype=float)
    n = len(y)
    out = []
    for pt in thresholds:
        positive = risk >= pt
        tp = ((positive == 1) & (y == 1)).sum()
        fp = ((positive == 1) & (y == 0)).sum()
        out.append(tp / n - fp / n * pt / (1 - pt))
    return np.asarray(out)


def decision_curve(d, cohort):
    rows = []
    for adjusted, label in [(False, "alpha3单独"), (True, "alpha3+年龄性别")]:
        X, y, state = _prepare(d, adjusted)
        if X is None:
            rows.append({"队列": cohort, "模型": label, "状态": state})
            continue
        try:
            model = _new_model().fit(X, y)
            risk = model.predict_proba(X)[:, 1]
            nb = _net_benefit(y, risk, THRESHOLDS)
            rng = np.random.default_rng(SEED + 20 + int(adjusted))
            boots = []
            for _ in range(BOOT_B):
                idx = rng.integers(0, len(y), len(y))
                if y.iloc[idx].nunique() < 2:
                    continue
                try:
                    bmodel = _new_model().fit(X.iloc[idx], y.iloc[idx])
                    brisk = bmodel.predict_proba(X.iloc[idx])[:, 1]
                    boots.append(_net_benefit(y.iloc[idx], brisk, THRESHOLDS))
                except Exception:
                    pass
            boot_array = np.asarray(boots)
            prevalence = y.mean()
            all_nb = prevalence - (1 - prevalence) * THRESHOLDS / (1 - THRESHOLDS)
            for j, pt in enumerate(THRESHOLDS):
                if len(boot_array) >= 100:
                    lo, hi = np.quantile(boot_array[:, j], [0.025, 0.975])
                else:
                    lo, hi = np.nan, np.nan
                rows.append({"队列": cohort, "模型": label, "n": len(y), "阈值概率": pt,
                             "净获益": nb[j], "净获益_bootstrap_CI低": lo,
                             "净获益_bootstrap_CI高": hi, "全部干预净获益": all_nb[j],
                             "无人干预净获益": 0.0, "相对全部干预净获益差": nb[j] - all_nb[j],
                             "有效bootstrap_n": len(boot_array),
                             "状态": "表观DCA及重拟合bootstrap；非外部验证"})
        except Exception as exc:
            rows.append({"队列": cohort, "模型": label, "状态": "DCA失败：%s" % type(exc).__name__})
    _save(pd.DataFrame(rows), "a38_IBS_DCA_净获益_bootstrap", cohort)


def alpha23_nnt(d, cohort):
    row = {"队列": cohort, "对照固定定义": "双低=alpha2<P67且alpha3<P67；双高=两者均>=P67；其余不进入该极端对照",
           "NNT定义": "1/(双低IBS风险-双高IBS风险)，仅描述性风险差倒数，不是治疗NNT"}
    needed = ["ibs_positive", "alpha2", "alpha3"]
    if any(c not in d for c in needed):
        row["状态"] = "所需变量缺失"
        _save(pd.DataFrame([row]), "a38_alpha2_alpha3_对照NNT", cohort)
        return
    sub = d[needed].dropna()
    if len(sub) < 20:
        row.update({"n": len(sub), "状态": "小样本，NNT不适用"})
        _save(pd.DataFrame([row]), "a38_alpha2_alpha3_对照NNT", cohort)
        return
    cut2, cut3 = sub["alpha2"].quantile(2 / 3), sub["alpha3"].quantile(2 / 3)
    low = ((sub["alpha2"] < cut2) & (sub["alpha3"] < cut3)).values
    high = ((sub["alpha2"] >= cut2) & (sub["alpha3"] >= cut3)).values
    y = sub["ibs_positive"].astype(int).values
    if low.sum() < 5 or high.sum() < 5:
        row.update({"n": len(sub), "双低n": int(low.sum()), "双高n": int(high.sum()),
                    "状态": "极端对照组小样本，NNT不适用"})
        _save(pd.DataFrame([row]), "a38_alpha2_alpha3_对照NNT", cohort)
        return
    risk_low, risk_high = y[low].mean(), y[high].mean()
    ard = risk_low - risk_high
    rng = np.random.default_rng(SEED + 30)
    diffs = []
    for _ in range(BOOT_B):
        il = rng.integers(0, low.sum(), low.sum())
        ih = rng.integers(0, high.sum(), high.sum())
        diffs.append(y[low][il].mean() - y[high][ih].mean())
    lo, hi = np.quantile(diffs, [0.025, 0.975])
    nnt = 1 / ard if ard > 0 else np.nan
    nnt_lo, nnt_hi = (1 / hi, 1 / lo) if lo > 0 and hi > 0 else (np.nan, np.nan)
    row.update({"n": len(sub), "alpha2_P67": cut2, "alpha3_P67": cut3,
                "双低n": int(low.sum()), "双高n": int(high.sum()),
                "双低IBS风险": risk_low, "双高IBS风险": risk_high,
                "绝对风险差_双低减双高": ard, "风险差_bootstrap_CI低": lo,
                "风险差_bootstrap_CI高": hi, "NNT_描述性": nnt,
                "NNT_CI低_仅风险差CI全正": nnt_lo, "NNT_CI高_仅风险差CI全正": nnt_hi,
                "状态": "已估计；风险差<=0或CI跨0时不报告有限NNT/区间"})
    _save(pd.DataFrame([row]), "a38_alpha2_alpha3_对照NNT", cohort)


def main():
    notes = pd.DataFrame([
        ["推断属性", "全部结果为探索性；除非另有可核查预注册，不得表述为确认性证据。"],
        ["PAR暴露", "固定为alpha3>=本队列P67（暴露）对其余；不根据结果切换暴露方向或切点。PAR为总体风险减非暴露风险后除以总体风险。"],
        ["DCA阈值", "预设0.01–0.50、步长0.01；净获益=TP/n-FP/n×pt/(1-pt)，同时列全部干预与无人干预。该范围不等于临床推荐阈值。"],
        ["bootstrap", "ROC摘要为表观评分重抽样区间；DCA在每个重抽样中重拟合模型；均不是外部验证。"],
        ["NNT边界", "α2+α3 NNT仅为固定双低/双高观察风险差倒数；不能解释为干预或治疗所需人数。"],
        ["敏感性排除", "精神用药=med_psych==1；主要共病固定为慢性痛、心脏病、高血压、糖尿病、高脂血症。"],
        ["不可识别", "任何亚组若IBS无病例-对照变异，明确标为不适用，不恢复或伪造ROC。"],
        ["隐私与资产", "仅输出聚合统计，不输出个人记录、题目级专有编码或映射。"]
    ], columns=["项目", "说明"])
    for cohort in C.COHORTS:
        d = C.load_cohort(cohort)
        if d is None:
            _save(pd.DataFrame([_status(cohort, "数据加载", "基底未就绪；未运行且不伪造结果")]),
                  "a38_运行状态", cohort)
            continue
        statuses = []
        jobs = [("PAR", lambda: population_attributable_fraction(d, cohort)),
                ("ROC/CV/bootstrap", lambda: roc_models(d, cohort)),
                ("用药共病敏感性", lambda: sensitivity(d, cohort)),
                ("DCA", lambda: decision_curve(d, cohort)),
                ("α2+α3对照NNT", lambda: alpha23_nnt(d, cohort))]
        for name, job in jobs:
            try:
                job()
                statuses.append(_status(cohort, name, "完成；具体不可估计项见各表状态"))
            except Exception as exc:
                statuses.append(_status(cohort, name, "跳过：%s：%s" % (type(exc).__name__, str(exc)[:160])))
        _save(pd.DataFrame(statuses), "a38_运行状态", cohort)
        _save(notes.assign(队列=cohort), "a38_分析口径与边界", cohort)
        print("[a38] %s IBS临床扩展完成（探索性）" % cohort)


if __name__ == "__main__":
    main()
