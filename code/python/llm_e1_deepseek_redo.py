# -*- coding: utf-8 -*-
"""E1 independent-reader comparison: public output-level downstream analysis.

Public inputs:
  code/data/derived/p3_base_{prolific,jianshu}.csv
  code/data/derived/p3_llm_predictions.csv

The prediction file contains only anonymous participant tokens and eight numeric
model estimates. Case text, item responses, raw model returns, response IDs,
input hashes, platform IDs, ID crosswalks, and API-call/build code remain in the
controlled research repository and are not part of this public package.

This script reproduces the five aggregate E1 tables from those public numeric
outputs. It does not rerun the model or reconstruct the proprietary backend.

Run:
  python code/python/llm_e1_deepseek_redo.py
"""
import os

import numpy as np
import pandas as pd

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
USER_ROOT = os.path.abspath(os.environ.get(
    "P3_ROOT", os.path.join(CODE_DIR, os.pardir, os.pardir)))
DERIVED = os.path.join(USER_ROOT, "code", "data", "derived")
OUT_DIR = os.path.abspath(os.environ.get(
    "P3_LLM_OUT", os.path.join(USER_ROOT, "data_output", "llm_e1")))
EXPECTED = {"prolific": 462, "jianshu": 504}
METRICS = ["alpha1", "alpha2", "alpha3", "alpha4", "coherence", "phq9", "gad7", "ibs"]


def _load_public(cohort):
    """Join anonymous model estimates to the public L2 base without any ID bridge."""
    base = pd.read_csv(
        os.path.join(DERIVED, "p3_base_%s.csv" % cohort),
        dtype=str, keep_default_na=False, encoding="utf-8-sig")
    pred = pd.read_csv(
        os.path.join(DERIVED, "p3_llm_predictions.csv"),
        dtype=str, keep_default_na=False, encoding="utf-8-sig")
    pred = pred[pred["cohort"] == cohort].copy()
    if len(base) != EXPECTED[cohort] or len(pred) != EXPECTED[cohort]:
        raise AssertionError("%s public LLM input dimensions are invalid" % cohort)
    if base["participant_token"].duplicated().any() or pred["participant_token"].duplicated().any():
        raise AssertionError("%s public LLM tokens are not unique" % cohort)

    direct = {
        "alpha1": "a_alpha1", "alpha2": "a_alpha2", "alpha3": "a_alpha3", "alpha4": "a_alpha4",
        "carrier_a": "a_f_a", "carrier_b": "a_f_b", "carrier_c": "a_f_c", "carrier_d": "a_f_d",
        "H_score": "a_H", "phq9_total": "phq9", "gad7_total": "gad7", "ibs_positive": "ibs",
        "chronic_pain": "chronic_pain", "maia_three_dim": "g_maia_three_dim",
        "aaq_total": "g_aaq_total", "dcs_total": "g_dcs_total", "bads_total": "g_bads_total",
    }
    out = base[["participant_token", "cohort", "data_source"] + list(direct)].copy()
    out = out.rename(columns=direct)
    out = out.merge(
        pred[["participant_token", "cohort"] + ["llm_" + m for m in METRICS]],
        on=["participant_token", "cohort"], how="inner", validate="one_to_one")
    if len(out) != EXPECTED[cohort]:
        raise AssertionError("%s public LLM prediction coverage is incomplete" % cohort)
    for col in list(direct.values()) + ["llm_" + m for m in METRICS]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    if out["a_H"].isna().any() or out[["llm_" + m for m in METRICS]].isna().any().any():
        raise AssertionError("%s public LLM required numeric values are missing" % cohort)
    return out


def _auc(y, s):
    y = np.asarray(y, float)
    s = np.asarray(s, float)
    m = ~(np.isnan(y) | np.isnan(s))
    y, s = y[m], s[m]
    if len(np.unique(y)) < 2:
        return np.nan
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(y, s)


def main():
    rows_corr, rows_auc, rows_gold, rows_dim, rows_combined = [], [], [], [], []
    for cohort in ["prolific", "jianshu"]:
        out = _load_public(cohort)
        print("[analyze/%s] merged public n=%d" % (cohort, len(out)))

        # Recovery correlations between independent-reader estimates and returned readouts.
        for factor in ["alpha1", "alpha2", "alpha3", "alpha4"]:
            pair = out[["llm_" + factor, "a_" + factor]].dropna()
            r = np.corrcoef(pair["llm_" + factor], pair["a_" + factor])[0, 1]
            rows_corr.append({
                "cohort": cohort, "pair": "llm_%s~algo_%s" % (factor, factor),
                "r": r, "n": len(pair)})
        pair = out[["llm_coherence", "a_H"]].dropna()
        r = np.corrcoef(pair["llm_coherence"], pair["a_H"])[0, 1]
        rows_corr.append({"cohort": cohort, "pair": "llm_coherence~algo_H", "r": r, "n": len(pair)})

        # Dimensionality of the four alpha estimates/readouts.
        def pc1(cols):
            x = out[cols].dropna().values
            x = (x - x.mean(0)) / x.std(0, ddof=1)
            ev = np.linalg.eigvals(np.corrcoef(x.T)).real
            return ev.max() / ev.sum()

        rows_dim.append({
            "cohort": cohort, "rating": "llm_alpha1-4",
            "PC1_pct": pc1(["llm_alpha1", "llm_alpha2", "llm_alpha3", "llm_alpha4"]) * 100})
        rows_dim.append({
            "cohort": cohort, "rating": "algo_alpha1-4",
            "PC1_pct": pc1(["a_alpha1", "a_alpha2", "a_alpha3", "a_alpha4"]) * 100})

        # Screening discrimination and exploratory same-sample combined models.
        for outcome, pred in [("phq9", "phq9"), ("gad7", "gad7"), ("ibs", "ibs")]:
            model_data = out[[outcome, "a_H", "llm_coherence", "llm_" + pred]].dropna().copy()
            y_med = (model_data[outcome] > model_data[outcome].median()).astype(int)
            y_cli = ((model_data[outcome] >= 10).astype(int) if outcome != "ibs"
                     else model_data[outcome].astype(int))
            for split, y in [("median", y_med), ("clinical" if outcome != "ibs" else "binary", y_cli)]:
                a_h = _auc(y, -model_data["a_H"])
                a_c = _auc(y, -model_data["llm_coherence"])
                a_phq = _auc(y, model_data["llm_" + pred])
                rows_auc.append({
                    "cohort": cohort, "outcome": outcome, "split": split,
                    "AUC_H": a_h, "AUC_llm_coherence": a_c,
                    "AUC_llm_est": a_phq, "n": int(y.notna().sum())})

                import statsmodels.api as sm
                for model_name, predictors in [
                        ("H+llm_coherence", ["a_H", "llm_coherence"]),
                        ("H+llm_direct", ["a_H", "llm_" + pred])]:
                    xx = sm.add_constant(model_data[predictors].astype(float), has_constant="add")
                    fit = sm.GLM(y.astype(int), xx, family=sm.families.Binomial()).fit()
                    if not bool(fit.converged):
                        raise RuntimeError("%s %s %s combined model did not converge" % (
                            cohort, outcome, model_name))
                    a_combined = _auc(y, fit.predict(xx))
                    rows_combined.append({
                        "cohort": cohort, "outcome": outcome, "split": split,
                        "model": model_name, "AUC_H": a_h, "AUC_combined": a_combined,
                        "delta_AUC_combined_minus_H": a_combined - a_h,
                        "n": int(y.notna().sum())})

        # External instruments are available only for the Prolific n=142 subsample.
        if cohort == "prolific":
            for col, name, direction in [
                    ("g_maia_three_dim", "MAIA-2", +1), ("g_aaq_total", "AAQ-II", -1),
                    ("g_dcs_total", "DCS", -1), ("g_bads_total", "BADS", +1)]:
                gold = out[["llm_coherence", "a_H", col]].dropna()
                r_llm = np.corrcoef(gold["llm_coherence"], gold[col])[0, 1]
                r_h = np.corrcoef(gold["a_H"], gold[col])[0, 1]
                rows_gold.append({
                    "criterion": name, "n": len(gold),
                    "r_llm_coherence": r_llm, "r_algo_H": r_h,
                    "theory_direction": direction,
                    "llm_dir_ok": int(np.sign(r_llm) == direction),
                    "algo_dir_ok": int(np.sign(r_h) == direction)})

    os.makedirs(OUT_DIR, exist_ok=True)
    pd.DataFrame(rows_corr).to_csv(
        os.path.join(OUT_DIR, "llm_e1_alpha_recovery.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(rows_auc).to_csv(
        os.path.join(OUT_DIR, "llm_e1_screening_auc.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(rows_gold).to_csv(
        os.path.join(OUT_DIR, "llm_e1_gold_standard.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(rows_dim).to_csv(
        os.path.join(OUT_DIR, "llm_e1_dimensionality.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(rows_combined).to_csv(
        os.path.join(OUT_DIR, "llm_e1_combined_auc.csv"), index=False, encoding="utf-8-sig")
    print("[analyze] five public aggregate tables written to %s" % OUT_DIR)


if __name__ == "__main__":
    main()
