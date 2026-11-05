# -*- coding: utf-8 -*-
"""a42_semantic_confound.py — 语义混淆矩阵、得分关联与结局文本核查。

内部模式：从用户指定目录的 BGE-M3 向量重算余弦矩阵，并与原矩阵逐格核验。
公开/轻量模式：向量不可用时只读取不含题目文本的相似度矩阵。

不在本模块重新编码文本，不输出题目全文、向量、逐人敏感数据或专有映射。
"""
import hashlib
import json
import os

import numpy as np
import pandas as pd
from scipy import stats

CODE_DIR = os.path.dirname(os.path.abspath(__file__))
USER_ROOT = os.path.abspath(os.environ.get("P3_ROOT", os.path.join(CODE_DIR, os.pardir, os.pardir)))
SOURCE_DIR = os.path.abspath(os.environ.get(
    "P3_SEMANTIC_DIR", os.path.join(USER_ROOT, "code", "data", "semantic")))
RESULT_PARENT = (os.path.join(os.path.abspath(os.environ["P3_RESULT_DIR"]), os.pardir)
                 if os.environ.get("P3_RESULT_DIR") else os.path.join(USER_ROOT, "data_output"))
OUT_DIR = os.path.abspath(os.environ.get(
    "P3_SEMANTIC_OUT", os.path.join(RESULT_PARENT, "a42_semantic_confound")))
GOLD_PREFIXES = {"MAIA": "alpha1", "AAQ": "alpha2", "DCS": "alpha3", "BADS": "alpha4"}


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _read_source_matrix(name):
    return pd.read_csv(os.path.join(SOURCE_DIR, name), index_col=0, encoding="utf-8-sig")


def _normalized(x):
    x = np.asarray(x, dtype=float)
    norm = np.linalg.norm(x, axis=1, keepdims=True)
    if not np.all(np.isfinite(x)) or np.any(norm <= 0):
        raise ValueError("嵌入矩阵含非有限值或零向量")
    return x / norm


def _stem_option_matrices():
    """优先从向量重算；无向量时读取公开安全的既有相似度矩阵。"""
    stem = _read_source_matrix("semantic_sim_stem.csv")
    option = _read_source_matrix("semantic_sim_options.csv")
    checks = {"embedding_recomputed": False}
    h_path = os.path.join(SOURCE_DIR, "hcf_embeddings.csv")
    g_path = os.path.join(SOURCE_DIR, "gold_embeddings.csv")
    o_path = os.path.join(SOURCE_DIR, "option_embeddings.csv")
    if all(os.path.exists(p) for p in (h_path, g_path, o_path)):
        h = _normalized(pd.read_csv(h_path, encoding="utf-8-sig").values)
        g = _normalized(pd.read_csv(g_path, encoding="utf-8-sig").values)
        o = _normalized(pd.read_csv(o_path, index_col=0, encoding="utf-8-sig").values)
        stem_calc = pd.DataFrame(h @ g.T, index=stem.index, columns=stem.columns)
        option_calc = pd.DataFrame(o @ g.T, index=option.index, columns=option.columns)
        np.testing.assert_allclose(stem_calc, stem, rtol=0, atol=5e-7)
        np.testing.assert_allclose(option_calc, option, rtol=0, atol=5e-7)
        checks.update({
            "embedding_recomputed": True,
            "stem_max_abs_diff": float(np.max(np.abs(stem_calc.values-stem.values))),
            "option_max_abs_diff": float(np.max(np.abs(option_calc.values-option.values))),
        })
        stem, option = stem_calc, option_calc
    if stem.shape != (20, 85) or option.shape != (100, 85):
        raise AssertionError("题干/选项矩阵形状错误: %s %s" % (stem.shape, option.shape))
    return stem, option, checks


def _summary_row(analysis, values, near=False):
    x = np.asarray(values, dtype=float).ravel()
    return {
        "analysis": analysis, "n": int(x.size), "mean": float(np.mean(x)),
        "min": float(np.min(x)), "max": float(np.max(x)),
        "n_gt_0_60": int(np.sum(x > .60)), "n_gt_0_70": int(np.sum(x > .70)),
        "n_gt_0_80": int(np.sum(x > .80)),
        "n_near_0_75_0_80": int(np.sum((x > .75) & (x <= .80))) if near else 0,
    }


def _correlation_row(analysis, x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    mask = np.isfinite(x) & np.isfinite(y)
    r, p = stats.pearsonr(x[mask], y[mask])
    return {"analysis": analysis, "r": r, "p": p, "n": int(mask.sum())}


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    stem, option, checks = _stem_option_matrices()
    outcome = _read_source_matrix("semantic_sim_phq_gad.csv")
    outcome_items = [c for c in outcome.columns if c.startswith("PHQ_") or c.startswith("GAD_")]
    if outcome.shape[0] != 20 or len(outcome_items) != 16:
        raise AssertionError("PHQ/GAD 相似度矩阵形状错误")

    # 公开安全的矩阵：仅题目/量表编号与相似度，不含文本或向量。
    stem.to_csv(os.path.join(OUT_DIR, "semantic_stem_matrix.csv"), encoding="utf-8-sig")
    option.to_csv(os.path.join(OUT_DIR, "semantic_option_matrix.csv"), encoding="utf-8-sig")

    summaries = [
        _summary_row("stem_all_20x85", stem.values),
        _summary_row("stem_per_item_max_over_85", stem.max(axis=1).values),
        _summary_row("option_all_100x85", option.values, near=True),
        _summary_row("option_per_option_max_over_85", option.max(axis=1).values),
        _summary_row("outcome_stem_all_20x16", outcome[outcome_items].values),
        _summary_row("outcome_stem_per_item_max_over_16", outcome[outcome_items].max(axis=1).values),
    ]
    pd.DataFrame(summaries).to_csv(os.path.join(OUT_DIR, "semantic_summary.csv"), index=False, encoding="utf-8-sig")

    scale_rows = []
    for level, matrix in (("stem", stem), ("option", option)):
        for prefix, alpha in GOLD_PREFIXES.items():
            cols = [c for c in matrix.columns if c.startswith(prefix + "_")]
            scale_rows.append({"level": level, "alpha": alpha, "gold_scale": prefix,
                               "n": int(matrix[cols].values.size),
                               "mean": float(matrix[cols].values.mean()),
                               "max": float(matrix[cols].values.max()),
                               "min": float(matrix[cols].values.min()),
                               "n_gt_0_80": int((matrix[cols].values > .80).sum())})
    pd.DataFrame(scale_rows).to_csv(os.path.join(OUT_DIR, "semantic_scale_summary.csv"), index=False, encoding="utf-8-sig")

    max_rows = []
    for level, matrix, id_name in (("stem", stem, "HCF_Item"), ("option", option, "Option")):
        for idx, row in matrix.iterrows():
            gold = row.idxmax()
            max_rows.append({"level": level, id_name: idx, "Max_Sim": row.max(), "Max_Gold_Item": gold})
    pd.DataFrame(max_rows).to_csv(os.path.join(OUT_DIR, "semantic_max_per_item.csv"), index=False, encoding="utf-8-sig")

    near = option.stack().rename("Semantic_Sim").reset_index()
    near.columns = ["Option", "Gold_Item", "Semantic_Sim"]
    near = near[(near["Semantic_Sim"] > .75) & (near["Semantic_Sim"] <= .80)].sort_values("Semantic_Sim", ascending=False)
    near.to_csv(os.path.join(OUT_DIR, "semantic_near_threshold_options.csv"), index=False, encoding="utf-8-sig")

    # 公开复现只读取聚合题目均分与选项选择率，不需要逐人作答。
    score_inputs = pd.read_csv(os.path.join(SOURCE_DIR, "semantic_score_inputs.csv"), encoding="utf-8-sig")
    if len(score_inputs) != 100 or score_inputs["HCF_Item"].nunique() != 20:
        raise AssertionError("语义聚合输入应为 20 题、100 选项")
    if not np.array_equal(score_inputs["Option"].astype(str).values, option.index.astype(str).values):
        raise AssertionError("语义聚合输入与选项矩阵顺序不一致")
    expected_items = np.repeat(stem.index.astype(str), 5)
    if not np.array_equal(score_inputs["HCF_Item"].astype(str).values, expected_items):
        raise AssertionError("语义聚合输入与题干矩阵的题目归属或顺序不一致")
    if (score_inputs.groupby("HCF_Item", sort=False)["item_mean_score"].nunique() != 1).any():
        raise AssertionError("同一题目的聚合均分不一致")
    item_means = score_inputs["item_mean_score"].to_numpy(dtype=float)[::5]
    endorsement = score_inputs["endorsement_rate"].to_numpy(dtype=float)
    if not np.isfinite(item_means).all() or not np.isfinite(endorsement).all():
        raise AssertionError("语义聚合输入含非有限值")
    corr_rows = [
        _correlation_row("stem_mean_similarity_vs_item_mean_score", item_means, stem.mean(axis=1).values),
        _correlation_row("stem_max_similarity_vs_item_mean_score", item_means, stem.max(axis=1).values),
        _correlation_row("option_mean_similarity_vs_endorsement_rate", endorsement, option.mean(axis=1).values),
        _correlation_row("option_max_similarity_vs_endorsement_rate", endorsement, option.max(axis=1).values),
    ]
    pd.DataFrame(corr_rows).to_csv(os.path.join(OUT_DIR, "semantic_score_correlations.csv"), index=False, encoding="utf-8-sig")
    # 公开复现用的聚合输入：仅题目均分与选项选择率，不含逐人作答。
    score_inputs = pd.DataFrame({
        "HCF_Item": np.repeat(np.asarray(stem.index), 5),
        "Option": np.asarray(option.index),
        "item_mean_score": np.repeat(item_means, 5),
        "endorsement_rate": np.asarray(endorsement, dtype=float),
    })
    score_inputs.to_csv(os.path.join(OUT_DIR, "semantic_score_inputs.csv"), index=False, encoding="utf-8-sig")

    manifest = {
        "model_recorded_in_source": "BAAI/bge-m3",
        "similarity": "L2-normalized cosine",
        "screening_threshold": 0.80,
        "near_threshold_rule": ">0.75 and <=0.80",
        "base_n": 462,
        "checks": checks,
        "source_files": {name: _sha256(os.path.join(SOURCE_DIR, name)) for name in (
            "semantic_sim_stem.csv", "semantic_sim_options.csv", "semantic_sim_phq_gad.csv")},
        "privacy": "No item text, embeddings, participant-level records, or proprietary mappings in outputs.",
    }
    with open(os.path.join(OUT_DIR, "semantic_validation.json"), "w", encoding="utf-8") as f:
        json.dump(manifest, f, ensure_ascii=False, indent=2)
    print("[a42] 输出目录: %s" % OUT_DIR)
    print(pd.DataFrame(summaries)[["analysis", "n", "mean", "min", "max", "n_gt_0_80"]].to_string(index=False))
    print(pd.DataFrame(corr_rows).to_string(index=False))


if __name__ == "__main__":
    main()
