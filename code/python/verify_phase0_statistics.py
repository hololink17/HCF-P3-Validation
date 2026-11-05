# verify_phase0_statistics.py —— C1/E3a/E2 独立复算与结果表核验
# 不导入被核验模块的统计函数；直接从派生基底和后台导出重算。
from pathlib import Path
import csv
import os

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm
from statsmodels.stats.nonparametric import rank_compare_2indep

ROOT = Path(os.environ.get("P3_ROOT", Path(__file__).resolve().parents[2])).resolve()
DERIVED = ROOT / "code" / "data" / "derived"
RESULT = ROOT / "data_output" / "analysis_results"
# 受控核验报告可写到仓库外，避免条目级受控结果覆盖或混入公开发布件。
OUT = (Path(os.environ["P0_REPORT_DIR"]).resolve() if os.environ.get("P0_REPORT_DIR")
       else RESULT / "P0_统计核验")
ROWS = []


def record(analysis, variable, check, observed, expected, tolerance=0.0):
    if isinstance(observed, (float, np.floating)) and isinstance(expected, (float, np.floating)):
        ok = bool(np.isfinite(observed) and abs(observed - expected) <= tolerance)
    else:
        ok = observed == expected
    ROWS.append({"analysis": analysis, "variable": variable, "check": check,
                 "observed": observed, "expected": expected, "tolerance": tolerance,
                 "status": "PASS" if ok else "FAIL"})
    if not ok:
        raise AssertionError("%s %s %s: observed=%r expected=%r" %
                             (analysis, variable, check, observed, expected))


def pooled_sd(a, b):
    return np.sqrt(((len(a)-1)*np.var(a, ddof=1)+(len(b)-1)*np.var(b, ddof=1)) /
                   (len(a)+len(b)-2))


def independent_tost(a, b, margin):
    diff = np.mean(a)-np.mean(b)
    se = np.sqrt(np.var(a, ddof=1)/len(a)+np.var(b, ddof=1)/len(b))
    df = (np.var(a, ddof=1)/len(a)+np.var(b, ddof=1)/len(b))**2 / (
        (np.var(a, ddof=1)/len(a))**2/(len(a)-1) +
        (np.var(b, ddof=1)/len(b))**2/(len(b)-1))
    p_lower = stats.t.sf((diff+margin)/se, df)
    p_upper = stats.t.cdf((diff-margin)/se, df)
    return diff, se, p_lower, p_upper, p_lower < .05 and p_upper < .05


def paired_tost(new, old, margin):
    d = np.asarray(new)-np.asarray(old)
    md = np.mean(d)
    se = np.std(d, ddof=1)/np.sqrt(len(d))
    p_lower = stats.t.sf((md+margin)/se, len(d)-1)
    p_upper = stats.t.cdf((md-margin)/se, len(d)-1)
    return md, p_lower, p_upper, p_lower < .05 and p_upper < .05


def verify_c1_e3a():
    long = pd.read_csv(DERIVED / "p3_base_phase2_long.csv")
    b2 = long.loc[(long.batch == 2) & (long.questionnaire == "新20")]
    reference = pd.read_csv(DERIVED / "p3_base_jianshu.csv")
    e3b_result = ROOT / "data_output" / "E3b_results"
    saved = pd.read_csv(e3b_result / "C1_预注册标准化主结果.csv").set_index("variable")
    for v in ["alpha1", "alpha2", "alpha3", "alpha4"]:
        a, b = b2[v].to_numpy(float), reference[v].to_numpy(float)
        fit = rank_compare_2indep(a, b, use_t=True)
        delta = 2*fit.prob1-1
        se = 2*np.sqrt(fit.var_prob)
        pl = stats.t.sf((delta+.147)/se, fit.df)
        pu = stats.t.cdf((delta-.147)/se, fit.df)
        eq = max(pl, pu) < .05
        record("C1", v, "n1+n2", len(a)+len(b), 744)
        record("C1", v, "cliff_delta", float(saved.loc[v, "cliff_delta"]), delta, 5.1e-7)
        record("C1", v, "se", float(saved.loc[v, "se"]), se, 5.1e-7)
        record("C1", v, "p_lower", float(saved.loc[v, "p_tost_lower"]), pl, 5.1e-7)
        record("C1", v, "p_upper", float(saved.loc[v, "p_tost_upper"]), pu, 5.1e-7)
        record("C1", v, "equivalent", bool(saved.loc[v, "equivalent"]), eq)
        record("C1", v, "注册判定", eq, v in ("alpha3", "alpha4"))

    wide = pd.read_csv(DERIVED / "p3_base_phase2_wide.csv")
    saved = pd.read_csv(e3b_result / "E3a_配对标准化敏感性.csv").set_index("variable")
    for v in ["alpha1", "alpha2", "alpha3", "alpha4"]:
        new, old = wide[v+"_new"].to_numpy(float), wide[v+"_old"].to_numpy(float)
        sd_diff = np.std(new-old, ddof=1)
        md, pl, pu, eq = paired_tost(new, old, .262*sd_diff)
        record("E3a", v, "n", len(new), 166)
        record("E3a", v, "raw_diff", float(saved.loc[v, "raw_diff"]), md, 5.1e-7)
        record("E3a", v, "sd_diff", float(saved.loc[v, "sd_diff"]), sd_diff, 5.1e-7)
        record("E3a", v, "standardized_margin", float(saved.loc[v, "standardized_margin"]),
               .262*sd_diff, 5.1e-7)
        record("E3a", v, "p_lower", float(saved.loc[v, "standardized_p_tost_lower"]), pl, 5.1e-7)
        record("E3a", v, "p_upper", float(saved.loc[v, "standardized_p_tost_upper"]), pu, 5.1e-7)
        record("E3a", v, "equivalent", bool(saved.loc[v, "standardized_equivalent"]), eq)
        record("E3a", v, "预定裁定", eq, v in ("alpha1", "alpha4"))


def read_dicts(path):
    for enc in ("utf-8-sig", "gb18030", "gbk"):
        try:
            with open(path, encoding=enc, newline="") as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            pass
    raise RuntimeError("无法解码：" + str(path))


def icc_stats(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    n, k = len(x), 2
    mat = np.column_stack([x, y])
    gm = mat.mean()
    ss_sub = k*np.sum((mat.mean(axis=1)-gm)**2)
    ss_for = n*np.sum((mat.mean(axis=0)-gm)**2)
    ss_res = np.sum((mat-gm)**2)-ss_sub-ss_for
    msr = ss_sub/(n-1)
    msf = ss_for/(k-1)
    mse = ss_res/((n-1)*(k-1))
    icc31 = (msr-mse)/(msr+(k-1)*mse)
    icc21 = (msr-mse)/(msr+(k-1)*mse+k*(msf-mse)/n)
    return icc31, icc21, stats.pearsonr(x, y).statistic, stats.spearmanr(x, y).statistic, y.mean()-x.mean()


def verify_e2():
    controlled_root = os.environ.get("P3_CONTROLLED_ROOT")
    if not controlled_root:
        raise RuntimeError("E2 受控核验需设置 P3_CONTROLLED_ROOT；公开 L2 不含原始配对桥接。")
    controlled_root = Path(controlled_root).resolve()
    exp1_path = controlled_root / "见数第一期" / "后台计算导出" / "见数第一期 第一套题 计算导出.csv"
    exp2_paths = [controlled_root / "见数第2期" / "后台计算导出" / f for f in
                  ["phase2_A40_先新后旧_旧题计算结果.csv", "phase2_A40_先旧后新_旧题计算结果.csv"]]
    raw1_paths = [controlled_root / "见数第一期" / "原始数据" / ("见数%d.csv" % i) for i in range(1, 8)]
    raw2_paths = [controlled_root / "见数第2期" / "原始数据" / f for f in
                  ["见数 第二期 40题 先第二套 66人.csv", "见数第二期40题 先第二套再第一套20人.csv",
                   "见数第二期40题先第一套后第二套30人.csv", "见数第二期40题先第一套后第二套50人.csv"]]

    def bridge(paths):
        out = {}
        for p in paths:
            for r in read_dicts(p):
                aid, uid = r["作答ID"].strip(), r["用户ID"].strip()
                if aid == "作答ID" and uid == "用户ID":  # 第二行Q编号，不是受访者
                    continue
                assert aid and uid and aid not in out
                out[aid] = uid
        return out

    b1, b2 = bridge(raw1_paths), bridge(raw2_paths)

    def score_map(rows, bridge_map):
        out = {}
        for r in rows:
            uid = bridge_map[r["external_profile_id"].strip()]
            assert uid not in out
            a = [float(r["scores_alpha%d" % i]) for i in range(1, 5)]
            h_field = os.environ.get("P3_H_FIELD", "h3_present")
            h = float(r[h_field])
            out[uid] = dict(zip(["alpha1", "alpha2", "alpha3", "alpha4", "H"], a+[h]))
        return out

    t1 = score_map(read_dicts(exp1_path), b1)
    rows2 = []
    for p in exp2_paths:
        rows2.extend(read_dicts(p))
    t2 = score_map(rows2, b2)
    common = sorted(set(t1) & set(t2))
    record("E2", "sample", "Time1_n", len(t1), 504)
    record("E2", "sample", "Time2_n", len(t2), 166)
    record("E2", "sample", "paired_n", len(common), 164)
    saved = pd.read_csv(ROOT / "data_output" / "E2_results" / "E2_icc_dual_20261031.csv").set_index("variable")
    alpha_reference = {"alpha1": .514, "alpha2": .571, "alpha3": .331, "alpha4": .600}
    for v in ["alpha1", "alpha2", "alpha3", "alpha4", "H"]:
        values = icc_stats([t1[u][v] for u in common], [t2[u][v] for u in common])
        for col, observed in zip(["ICC_3_1", "ICC_2_1", "Pearson_r", "Spearman_rho", "mean_T2_minus_T1"], values):
            tolerance = .00051 if col != "mean_T2_minus_T1" else .000051
            record("E2", v, col, float(saved.loc[v, col]), float(observed), tolerance)
        if v in alpha_reference:
            record("E2", v, "α原裁定值不变", float(saved.loc[v, "ICC_3_1"]),
                   alpha_reference[v], .00051)


def _matrix(d, cols):
    return np.column_stack([np.ones(len(d)), d[cols].to_numpy(dtype=float)])


def _ols_fit(y, x):
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    resid = y - x @ beta
    sse = float(resid @ resid)
    centered = y - y.mean()
    r2 = float(1 - sse / (centered @ centered))
    df_resid = len(y) - np.linalg.matrix_rank(x)
    return {"resid": resid, "sse": sse, "r2": r2, "df_resid": df_resid}


def _compare_nested(d, outcome, base_cols, added_cols):
    y = d[outcome].to_numpy(dtype=float)
    f0 = _ols_fit(y, _matrix(d, base_cols))
    f1 = _ols_fit(y, _matrix(d, base_cols + added_cols))
    df_num = f0["df_resid"] - f1["df_resid"]
    f = ((f0["sse"] - f1["sse"]) / df_num) / (f1["sse"] / f1["df_resid"])
    p = float(stats.f.sf(f, df_num, f1["df_resid"]))
    return f0, f1, float(f), df_num, p


def _delong_difference(y, score0, score1):
    """独立结构分量复算两相关ROC的DeLong差值检验。"""
    pos, neg = y == 1, y == 0
    m, n = int(pos.sum()), int(neg.sum())
    components10, components01, aucs = [], [], []
    for score in (score0, score1):
        xp, xn = score[pos], score[neg]
        comparisons = ((xp[:, None] > xn[None, :]).astype(float) +
                       .5 * (xp[:, None] == xn[None, :]))
        components10.append(comparisons.mean(axis=1))
        components01.append(comparisons.mean(axis=0))
        aucs.append(float(comparisons.mean()))
    cov10 = np.cov(np.vstack(components10), ddof=1)
    cov01 = np.cov(np.vstack(components01), ddof=1)
    variance = float(cov10[0, 0] / m + cov10[1, 1] / m - 2 * cov10[0, 1] / m +
                     cov01[0, 0] / n + cov01[1, 1] / n - 2 * cov01[0, 1] / n)
    delta = aucs[1] - aucs[0]
    z = delta / np.sqrt(variance)
    return aucs, delta, float(z), float(2 * stats.norm.sf(abs(z)))


def _profile_permutation(d, outcome, base_cols, added_cols, observed_delta, seed=20261101, B=2000):
    y = d[outcome].to_numpy(dtype=float)
    x0, x1 = _matrix(d, base_cols), _matrix(d, base_cols + added_cols)
    fit0 = _ols_fit(y, x0)
    fitted0 = y - fit0["resid"]

    def relative_delta(yy):
        r0 = _ols_fit(yy, x0)["resid"]
        r1 = _ols_fit(yy, x1)["resid"]
        return float(1 - (r1 @ r1) / (r0 @ r0))

    observed = relative_delta(y)
    assert np.isclose(observed, observed_delta / (1 - fit0["r2"]), atol=1e-12, rtol=1e-10)
    rng = np.random.default_rng(seed)
    count = 0
    for _ in range(B):
        count += relative_delta(fitted0 + rng.permutation(fit0["resid"])) >= observed - 1e-12
    return (count + 1) / (B + 1)


def verify_a41():
    """不导入a41/common；由Gold基底独立复算五张聚合结果表。"""
    g = pd.read_csv(DERIVED / "p3_base_gold.csv")
    outdir = RESULT / "a41_gold_incremental_ipsative" / "prolific"
    alphas = ["alpha1", "alpha2", "alpha3", "alpha4"]
    golds = ["maia_three_dim", "aaq_total", "dcs_total", "bads_total"]
    symptoms = ["phq9_total", "gad7_total"]
    record("a41", "Gold", "n", len(g), 142)

    saved = pd.read_csv(outdir / "a41_Gold增量效度.csv").set_index("outcome")
    for outcome in symptoms:
        f0, f1, f, df_num, p = _compare_nested(g, outcome, golds, alphas)
        row = saved.loc[outcome]
        for col, expected in [("R2_base", f0["r2"]), ("R2_full", f1["r2"]),
                              ("delta_R2", f1["r2"] - f0["r2"]), ("F", f),
                              ("df_numerator", df_num), ("p", p)]:
            record("a41增量", outcome, col, float(row[col]), expected, 1e-10)

        y = (g[outcome] >= 10).astype(int).to_numpy()
        fits = [sm.GLM(y, _matrix(g, cols), family=sm.families.Binomial()).fit()
                for cols in (golds, golds + alphas)]
        record("a41临床增量", outcome, "两模型均收敛",
               all(bool(fit.converged) for fit in fits), True)
        aucs, delta_auc, z, p_delong = _delong_difference(
            y, fits[0].predict(), fits[1].predict())
        row = saved.loc[outcome + "_clinical10"]
        pseudo = [1 - fit.deviance / fit.null_deviance for fit in fits]
        for col, expected in [("R2_base", pseudo[0]), ("R2_full", pseudo[1]),
                              ("delta_R2", pseudo[1] - pseudo[0]),
                              ("AUC_base", aucs[0]), ("AUC_full", aucs[1]),
                              ("delta_AUC", delta_auc), ("DeLong_z", z), ("p", p_delong)]:
            record("a41临床增量", outcome, col, float(row[col]), expected, 1e-10)

    saved_partial = pd.read_csv(outdir / "a41_症状偏相关_16对.csv")
    for _, row in saved_partial.iterrows():
        x = _matrix(g, symptoms)
        ra = _ols_fit(g[row["alpha"]].to_numpy(float), x)["resid"]
        rg = _ols_fit(g[row["criterion"]].to_numpy(float), x)["resid"]
        r = float(np.corrcoef(ra, rg)[0, 1])
        df = len(g) - len(symptoms) - 2
        t = r * np.sqrt(df / (1 - r * r))
        p = float(2 * stats.t.sf(abs(t), df))
        key = row["criterion"] + "/" + row["alpha"]
        record("a41偏相关", key, "partial_r", float(row["partial_r"]), r, 1e-12)
        record("a41偏相关", key, "df", int(row["df"]), df)
        record("a41偏相关", key, "p", float(row["p"]), p, 1e-12)
        record("a41偏相关", key, "Bonferroni", float(row["p_bonferroni_16"]),
               min(1, 16 * p), 1e-12)

    saved_joint = pd.read_csv(outdir / "a41_组内剖面联合检验.csv").set_index(["outcome", "model"])
    saved_coef = pd.read_csv(outdir / "a41_组内剖面系数.csv")
    saved_corr = pd.read_csv(outdir / "a41_组内剖面单相关.csv")
    record("a41剖面", "joint_table", "行数", len(saved_joint), 8)
    record("a41剖面", "coefficient_table", "行数", len(saved_coef), 12)
    record("a41剖面", "correlation_table", "行数", len(saved_corr), 16)
    for outcome in golds:
        centered4 = g[alphas].sub(g[alphas].mean(axis=1), axis=0)
        centered4.columns = ["profile_" + c for c in alphas]
        work4 = pd.concat([g[[outcome] + symptoms], centered4], axis=1)
        added4 = list(centered4.columns[:-1])
        f0, f1, f, df_num, p = _compare_nested(work4, outcome, symptoms, added4)
        row = saved_joint.loc[(outcome, "full_4alpha_profile")]
        for col, expected in [("R2_base", f0["r2"]), ("R2_full", f1["r2"]),
                              ("delta_R2", f1["r2"] - f0["r2"]), ("F", f),
                              ("df_numerator", df_num), ("p", p)]:
            record("a41剖面", outcome + "/full", col, float(row[col]), expected, 1e-10)
        record("a41剖面", outcome + "/full", "Bonferroni",
               float(row["p_bonferroni_4criteria"]), min(1, 4 * p), 1e-12)
        record("a41剖面", outcome + "/full", "置换方法", row["permutation_method"], "Freedman-Lane")
        perm = _profile_permutation(work4, outcome, symptoms, added4, f1["r2"] - f0["r2"])
        record("a41剖面", outcome + "/full", "置换p", float(row["permutation_p"]), perm, 1e-15)

        coef_sub = work4[[outcome] + symptoms + added4].dropna()
        coef_design = sm.add_constant(coef_sub[symptoms + added4].astype(float),
                                      has_constant="add")
        coef_fit = sm.OLS(coef_sub[outcome].astype(float), coef_design).fit()
        coef_rows = saved_coef.loc[saved_coef["criterion"] == outcome].set_index("term")
        record("a41剖面系数", outcome, "模型标签唯一",
               set(saved_coef.loc[saved_coef["criterion"] == outcome, "model"]),
               {"full_4alpha_profile"})
        for term in added4:
            coef_row = coef_rows.loc[term]
            key = outcome + "/" + term
            for col, expected in [("beta", coef_fit.params[term]), ("SE", coef_fit.bse[term]),
                                  ("t", coef_fit.tvalues[term]), ("p", coef_fit.pvalues[term])]:
                record("a41剖面系数", key, col, float(coef_row[col]), float(expected), 1e-12)
            record("a41剖面系数", key, "n", int(coef_row["n"]), len(coef_sub))

        corr_sub = work4[[outcome] + symptoms + list(centered4.columns)].dropna()
        corr_rows = saved_corr.loc[saved_corr["criterion"] == outcome].set_index("term")
        for term in centered4.columns:
            corr_row = corr_rows.loc[term]
            expected_r, expected_p = stats.pearsonr(corr_sub[term], corr_sub[outcome])
            key = outcome + "/" + term
            record("a41剖面单相关", key, "r", float(corr_row["r"]), float(expected_r), 1e-12)
            record("a41剖面单相关", key, "p", float(corr_row["p"]), float(expected_p), 1e-12)
            record("a41剖面单相关", key, "Bonferroni",
                   float(corr_row["p_bonferroni_16"]), min(1, 16 * float(expected_p)), 1e-12)
            record("a41剖面单相关", key, "n", int(corr_row["n"]), len(corr_sub))

        cols3 = ["alpha1", "alpha2", "alpha4"]
        centered3 = g[cols3].sub(g[cols3].mean(axis=1), axis=0)
        centered3.columns = ["profile_" + c for c in cols3]
        work3 = pd.concat([g[[outcome] + symptoms], centered3], axis=1)
        added3 = list(centered3.columns[:-1])
        f0, f1, f, df_num, p = _compare_nested(work3, outcome, symptoms, added3)
        row = saved_joint.loc[(outcome, "exclude_alpha3_recentered")]
        for col, expected in [("R2_base", f0["r2"]), ("R2_full", f1["r2"]),
                              ("delta_R2", f1["r2"] - f0["r2"]), ("F", f),
                              ("df_numerator", df_num), ("p", p)]:
            record("a41剖面", outcome + "/exclude3", col, float(row[col]), expected, 1e-10)
        record("a41剖面", outcome + "/exclude3", "Bonferroni",
               float(row["p_bonferroni_4criteria"]), min(1, 4 * p), 1e-12)
        perm = _profile_permutation(work3, outcome, symptoms, added3, f1["r2"] - f0["r2"])
        record("a41剖面", outcome + "/exclude3", "置换p", float(row["permutation_p"]), perm, 1e-15)


def _auc_ci(y, score):
    """由阳性/阴性比较矩阵独立计算AUC及DeLong正态CI。"""
    y, score = np.asarray(y, float), np.asarray(score, float)
    pos, neg = score[y == 1], score[y == 0]
    comparisons = ((pos[:, None] > neg[None, :]).astype(float) +
                   .5 * (pos[:, None] == neg[None, :]))
    v10, v01 = comparisons.mean(axis=1), comparisons.mean(axis=0)
    auc = float(comparisons.mean())
    variance = np.var(v10, ddof=1)/len(pos) + np.var(v01, ddof=1)/len(neg)
    half = stats.norm.ppf(.975)*np.sqrt(variance)
    return auc, max(0.0, auc-half), min(1.0, auc+half)


def _bh(values):
    """独立实现BH单调化，不调用被核验模块或multipletests。"""
    values = np.asarray(values, float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked*len(values)/np.arange(1, len(values)+1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty(len(values), float)
    out[order] = np.minimum(adjusted, 1.0)
    return out


def verify_a04():
    """双队列144个单变量、144个校正及126个症状增量模型逐格独立复算。"""
    alphas = ["alpha1", "alpha2", "alpha3", "alpha4"]
    covars = ["age_group", "gender", "education"]
    symptoms = ["phq9_total", "gad7_total"]
    for cohort in ["prolific", "jianshu"]:
        d = pd.read_csv(DERIVED / ("p3_base_%s.csv" % cohort))
        outdir = RESULT / "a04_predictive" / cohort

        saved = pd.read_csv(outdir / "a04_九读数疾病perSD.csv")
        record("a04单变量", cohort, "模型数", len(saved), 72)
        record("a04单变量", cohort, "全部收敛", bool((saved["状态"] == "ok").all()), True)
        for _, row in saved.iterrows():
            outcome, pred = row["结局"], row["预测因子"]
            key = "%s/%s/%s" % (cohort, outcome, pred)
            sub = d[[outcome, pred]].dropna().copy()
            y = sub[outcome].astype(float)
            z = (sub[pred]-sub[pred].mean())/sub[pred].std(ddof=1)
            x = sm.add_constant(z.astype(float), has_constant="add")
            fit = sm.GLM(y, x, family=sm.families.Binomial()).fit()
            record("a04单变量", key, "收敛", bool(fit.converged), True)
            ci = fit.conf_int().loc[pred]
            auc = _auc_ci(y, -sub[pred].to_numpy(float))
            expected = {"AUC": auc[0], "AUC_CI低": auc[1], "AUC_CI高": auc[2],
                        "OR_per_SD": np.exp(fit.params[pred]),
                        "OR_CI低": np.exp(ci[0]), "OR_CI高": np.exp(ci[1]),
                        "p": fit.pvalues[pred]}
            for col, value in expected.items():
                record("a04单变量", key, col, float(row[col]), float(value), 1e-10)
            record("a04单变量", key, "n", int(row["n"]), len(sub))
            record("a04单变量", key, "阳性数", int(row["阳性数"]), int(y.sum()))

        saved = pd.read_csv(outdir / "a04_疾病OR_校正.csv")
        record("a04校正", cohort, "模型数", len(saved), 72)
        record("a04校正", cohort, "全部收敛", bool((saved["状态"] == "ok").all()), True)
        independent_p = []
        for _, row in saved.iterrows():
            outcome, pred = row["结局"], row["预测因子"]
            key = "%s/%s/%s" % (cohort, outcome, pred)
            sub = d[[outcome, pred]+covars].dropna().copy()
            x = pd.get_dummies(sub[[pred]+covars], columns=covars,
                               drop_first=True, dtype=float)
            x = sm.add_constant(x.astype(float), has_constant="add")
            fit = sm.GLM(sub[outcome].astype(float), x,
                         family=sm.families.Binomial()).fit()
            record("a04校正", key, "收敛", bool(fit.converged), True)
            ci = fit.conf_int().loc[pred]
            independent_p.append(float(fit.pvalues[pred]))
            expected = {"校正OR": np.exp(fit.params[pred]),
                        "CI低": np.exp(ci[0]), "CI高": np.exp(ci[1]),
                        "校正OR_per_0_1": np.exp(.1*fit.params[pred]),
                        "CI低_per_0_1": np.exp(.1*ci[0]),
                        "CI高_per_0_1": np.exp(.1*ci[1]), "p": fit.pvalues[pred]}
            for col, value in expected.items():
                record("a04校正", key, col, float(row[col]), float(value), 1e-10)
            record("a04校正", key, "n", int(row["n"]), len(sub))
        expected_q9 = _bh(independent_p)
        alpha_mask = saved["预测因子"].isin(alphas).to_numpy()
        expected_q4 = np.full(len(saved), np.nan)
        expected_q4[alpha_mask] = _bh(np.asarray(independent_p)[alpha_mask])
        for i, row in saved.iterrows():
            key = "%s/%s/%s" % (cohort, row["结局"], row["预测因子"])
            record("a04校正FDR", key, "q_BH_9readouts_8outcomes",
                   float(row["q_BH_9readouts_8outcomes"]), float(expected_q9[i]), 1e-12)
            if alpha_mask[i]:
                record("a04校正FDR", key, "q_BH_4alpha_8outcomes",
                       float(row["q_BH_4alpha_8outcomes"]), float(expected_q4[i]), 1e-12)

        saved = pd.read_csv(outdir / "a04_控症状后疾病增量.csv")
        record("a04症状增量", cohort, "模型数", len(saved), 63)
        record("a04症状增量", cohort, "全部收敛", bool((saved["状态"] == "ok").all()), True)
        independent_lr = []
        for _, row in saved.iterrows():
            outcome, pred = row["结局"], row["预测因子"]
            key = "%s/%s/%s" % (cohort, outcome, pred)
            sub = d[[outcome, pred]+covars+symptoms].dropna().copy()
            y = sub[outcome].astype(float)
            x = pd.get_dummies(sub[covars], columns=covars, drop_first=True, dtype=float)
            for symptom in symptoms:
                x[symptom] = sub[symptom].astype(float)
            x0 = sm.add_constant(x.astype(float), has_constant="add")
            x1 = x0.copy()
            x1[pred] = sub[pred].astype(float)
            m0 = sm.GLM(y, x0, family=sm.families.Binomial()).fit()
            m1 = sm.GLM(y, x1, family=sm.families.Binomial()).fit()
            record("a04症状增量", key, "两模型均收敛",
                   bool(m0.converged and m1.converged), True)
            lr = max(0.0, 2*(m1.llf-m0.llf))
            lr_p = float(stats.chi2.sf(lr, 1))
            independent_lr.append(lr_p)
            auc0, auc1 = _auc_ci(y, m0.predict())[0], _auc_ci(y, m1.predict())[0]
            ci = m1.conf_int().loc[pred]
            expected = {"AUC_base": auc0, "AUC_full": auc1, "delta_AUC": auc1-auc0,
                        "LR_chi2": lr, "LR_p": lr_p, "coef_p": m1.pvalues[pred],
                        "OR_per_0_1": np.exp(.1*m1.params[pred]),
                        "CI低_per_0_1": np.exp(.1*ci[0]),
                        "CI高_per_0_1": np.exp(.1*ci[1])}
            for col, value in expected.items():
                record("a04症状增量", key, col, float(row[col]), float(value), 1e-10)
            record("a04症状增量", key, "n", int(row["n"]), len(sub))
        saved["_独立LR_p"] = independent_lr
        for outcome, indices in saved.groupby("结局").groups.items():
            idx = list(indices)
            expected_q = _bh(saved.loc[idx, "_独立LR_p"])
            for i, q in zip(idx, expected_q):
                row = saved.loc[i]
                key = "%s/%s/%s" % (cohort, outcome, row["预测因子"])
                record("a04症状FDR", key, "q_BH_within_outcome_9",
                       float(row["q_BH_within_outcome_9"]), float(q), 1e-12)


def verify_a25():
    covars = ["age_group", "gender", "education"]
    for cohort in ["prolific", "jianshu"]:
        d = pd.read_csv(DERIVED / ("p3_base_%s.csv" % cohort))
        saved = pd.read_csv(RESULT / "a25_replication_extension" / cohort /
                            "a25_BMI四alpha增量回归.csv").set_index("alpha")
        for alpha in ["alpha1", "alpha2", "alpha3", "alpha4"]:
            key = "%s/%s" % (cohort, alpha)
            sub = d[["bmi", alpha]+covars].dropna().copy()
            x0 = pd.get_dummies(sub[covars], columns=covars, drop_first=True, dtype=float)
            x0 = sm.add_constant(x0, has_constant="add")
            x1 = x0.copy()
            x1[alpha] = sub[alpha].astype(float)
            m0 = sm.OLS(sub["bmi"].astype(float), x0).fit()
            m1 = sm.OLS(sub["bmi"].astype(float), x1).fit()
            row = saved.loc[alpha]
            expected = {"beta": m1.params[alpha], "SE": m1.bse[alpha],
                        "t": m1.tvalues[alpha], "p": m1.pvalues[alpha],
                        "delta_R2": m1.rsquared-m0.rsquared, "base_R2": m0.rsquared}
            for col, value in expected.items():
                record("a25_BMI", key, col, float(row[col]), float(value), 1e-10)
            record("a25_BMI", key, "n", int(row["n"]), len(sub))


def verify_a27():
    controlled_root = os.environ.get("P3_CONTROLLED_ROOT")
    if not controlled_root:
        raise RuntimeError("a27 条目级受控核验需设置 P3_CONTROLLED_ROOT；公开 L2 不含 PHQ-9 逐题字段。")
    controlled_derived = Path(controlled_root).resolve() / "代码" / "data" / "derived"
    for cohort in ["prolific", "jianshu", "gold"]:
        d = pd.read_csv(controlled_derived / ("p3_base_%s.csv" % cohort)).copy()
        d["si"] = np.where(d["phq_44"].notna(), (d["phq_44"] >= 1).astype(float), np.nan)
        other8_cols = ["phq_%d" % i for i in range(36, 44)]
        d["phq_other8"] = d[other8_cols].sum(axis=1, min_count=8)
        outdir = Path(controlled_root).resolve() / "分析结果" / "a27_suicidal_ideation" / cohort
        saved = pd.read_csv(outdir / "a27_条目描述与四alpha.csv").set_index("alpha")
        for alpha in ["alpha1", "alpha2", "alpha3", "alpha4"]:
            key = "%s/%s" % (cohort, alpha)
            sub = pd.DataFrame({"item": d["phq_44"], "si": d["si"],
                                "alpha": d[alpha], "other8": d["phq_other8"]}).dropna()
            pos = sub.loc[sub.si == 1, "alpha"].to_numpy(float)
            neg = sub.loc[sub.si == 0, "alpha"].to_numpy(float)
            pooled = np.sqrt(((len(pos)-1)*pos.var(ddof=1)+(len(neg)-1)*neg.var(ddof=1)) /
                             (len(pos)+len(neg)-2))
            rho, rho_p = stats.spearmanr(sub["alpha"], sub["item"])
            x = _matrix(sub, ["other8"])
            ra = _ols_fit(sub["alpha"].to_numpy(float), x)["resid"]
            ri = _ols_fit(sub["item"].to_numpy(float), x)["resid"]
            r = float(np.corrcoef(ra, ri)[0, 1])
            df = len(sub)-3
            partial_p = float(2*stats.t.sf(abs(r*np.sqrt(df/(1-r*r))), df))
            row = saved.loc[alpha]
            expected = {"SI阳性数": int(sub.si.sum()), "SI阳性率": sub.si.mean(),
                        "阳性组均值": pos.mean(), "阴性组均值": neg.mean(),
                        "均值差_阳性减阴性": pos.mean()-neg.mean(),
                        "Cohen_d": (pos.mean()-neg.mean())/pooled,
                        "Welch_p": stats.ttest_ind(pos, neg, equal_var=False).pvalue,
                        "MannWhitney_p": stats.mannwhitneyu(pos, neg, alternative="two-sided").pvalue,
                        "Spearman_r": rho, "Spearman_p": rho_p,
                        "partial_r_控制其余8题": r, "partial_df": df, "partial_p": partial_p}
            for col, value in expected.items():
                observed = row[col]
                if isinstance(value, (float, np.floating)):
                    record("a27条目", key, col, float(observed), float(value), 1e-10)
                else:
                    record("a27条目", key, col, int(observed), int(value))
            record("a27条目", key, "n", int(row["n"]), len(sub))

        saved_inc = pd.read_csv(outdir / "a27_alpha1增量效度.csv")
        for _, row in saved_inc.iterrows():
            base = ["phq_other8"] if row["模型"] == "PHQ其余8题+alpha1" else ["maia_three_dim"]
            key = "%s/%s" % (cohort, row["模型"])
            dd = d[["si"]+base+["alpha1"]].dropna()
            x0 = sm.add_constant(dd[base].astype(float), has_constant="add")
            x1 = sm.add_constant(dd[base+["alpha1"]].astype(float), has_constant="add")
            m0 = sm.GLM(dd.si.astype(float), x0, family=sm.families.Binomial()).fit()
            m1 = sm.GLM(dd.si.astype(float), x1, family=sm.families.Binomial()).fit()
            record("a27增量", key, "两模型均收敛", bool(m0.converged and m1.converged), True)
            lr = max(0.0, 2*(m1.llf-m0.llf))
            expected = {"alpha1_logOR": m1.params["alpha1"],
                        "alpha1_OR": np.exp(m1.params["alpha1"]),
                        "alpha1_p": m1.pvalues["alpha1"],
                        "增量LR_p": stats.chi2.sf(lr, 1)}
            for col, value in expected.items():
                record("a27增量", key, col, float(row[col]), float(value), 1e-10)
            record("a27增量", key, "n", int(row["n"]), len(dd))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    report_path = OUT / "P0_独立核验.csv"
    public_l2 = not bool(os.environ.get("P3_CONTROLLED_ROOT"))
    steps = [
        ("C1_E3a", verify_c1_e3a, False),
        ("E2", verify_e2, True),
        ("a41", verify_a41, False),
        ("a04", verify_a04, False),
        ("a25", verify_a25, False),
        ("a27", verify_a27, True),
    ]
    for name, function, controlled_only in steps:
        if public_l2 and controlled_only:
            ROWS.append({"analysis": name, "variable": "controlled_only",
                         "check": "公开L2执行边界", "observed": "SKIPPED",
                         "expected": "受控L3原始桥接或逐题数据", "tolerance": "",
                         "status": "SKIPPED_CONTROLLED_L3_REQUIRED"})
            print("[P0] %s 需要受控L3数据，公开L2模式明确跳过。" % name)
            continue
        try:
            function()
        except Exception as exc:
            ROWS.append({"analysis": name, "variable": "runner", "check": "模块执行",
                         "observed": "%s: %s" % (type(exc).__name__, exc),
                         "expected": "无异常", "tolerance": "", "status": "FAIL"})
            pd.DataFrame(ROWS).to_csv(report_path, index=False, encoding="utf-8-sig")
            print("P0核验在 %s 中断；失败与已完成检查已写入：%s" % (name, report_path))
            raise
    report = pd.DataFrame(ROWS)
    report.to_csv(report_path, index=False, encoding="utf-8-sig")
    passed = int((report["status"] == "PASS").sum())
    skipped = int((report["status"] == "SKIPPED_CONTROLLED_L3_REQUIRED").sum())
    print("P0核验完成：%d项PASS，%d项受控跳过，0项FAIL；报告：%s" % (passed, skipped, report_path))
    if public_l2:
        print("公开L2模式未执行E2原始桥接及a27逐题重算；不得据此声称受控血缘或条目级核验已完成。")
    else:
        print("受控模式已执行C1/E3a、E2、a41五表、a04、a25及a27的独立复算。")


if __name__ == "__main__":
    main()
