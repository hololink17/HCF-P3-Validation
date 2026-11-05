# =============================================================
# common.py —— P3 工具验证 · Python 移植单一真值层
# -------------------------------------------------------------
# 移植自 R 管线：
#   · 00_shared_rules.R      （跨平台共用规则，唯一真值）
#   · 10_analysis_00_common.R（分析层公共入口）
#
# 本文件是「唯一口径真值」，所有 aNN_*.py 模块都 import 本文件，
# 禁止各模块自行定义计分/分型/派生规则，防止口径分叉。
#
# 数据来源：code/data/derived/p3_base_{prolific,jianshu,gold}.csv。
# 公开分析只读取后台直接输出的 α1–α4、carrier_a–d 与 H_score，
# 以及外部变量；不从原始作答重构或复算任何内部分数。
# 唯一额外 HCF 派生量为书籍已公开的 H_trend = alpha1^2，仅属下游变换。
# 外部公开量表按各自公开计分规则处理。
#
# 版本标识：问卷 科研-第一套题 ｜ 算法 算法v1.1。
# =============================================================

import os
import numpy as np
import pandas as pd

# ---------- 路径（唯一真值；默认相对本文件，可用环境变量覆盖）----------
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
USER_ROOT = os.path.abspath(os.environ.get(
    "P3_ROOT", os.path.join(CODE_DIR, os.pardir, os.pardir)))
PROJ = USER_ROOT  # 公开包根目录；可用 P3_ROOT 指向受控完整数据根目录
DERIVED = os.path.join(USER_ROOT, "code", "data", "derived")
RESULT_DIR = os.path.abspath(os.environ.get(
    "P3_RESULT_DIR", os.path.join(USER_ROOT, "data_output", "analysis_results")))

PROLIFIC_F = os.path.join(DERIVED, "p3_base_prolific.csv")
JIANSHU_F = os.path.join(DERIVED, "p3_base_jianshu.csv")
GOLD_F = os.path.join(DERIVED, "p3_base_gold.csv")
PROLIFIC_RAW_CLEAN = os.path.join(USER_ROOT, "Prolific", "原始数据", "prolific_clean_人口学PHQ9疾病.csv")

Q_VERSION = "科研-第一套题"
ALGO_VERSION = "算法v1.1"
VERSION_TAG = "问卷与映射 %s ｜ 算法参数 %s" % (Q_VERSION, ALGO_VERSION)

COHORTS = ["prolific", "jianshu"]
ALPHA_AC = ["alpha1", "alpha2", "alpha3", "alpha4"]

# 数值化列（与 R numc 对齐 + 题目级）
_NUM_COLS = [
    "alpha1", "alpha2", "alpha3", "alpha4",
    "carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score",
    "phq9_total", "gad7_total", "gad7_item_count",
    "ibs_positive", "ibs_symptom_count",
    "chronic_pain", "heart_disease", "hypertension", "diabetes", "hyperlipidemia",
    "bmi", "obese", "height_cm", "weight_kg", "age", "sleep_hour", "smoke", "alcohol_freq", "alcohol_binary",
    "med_hypertension", "med_diabetes", "med_lipid", "med_psych",
    "sleep_disorder_positive", "sleep_symptom_count",
    "anchor1", "anchor2", "anchor3", "anchor4",
    "maia_noticing", "maia_not_distracting", "maia_not_worrying", "maia_three_dim",
    "aaq_total", "dcs_total", "bads_total",
    "duration_sec", "depression", "phq9_normal",
    "attention_pass", "inattentive",
] + ["hcf_%d" % i for i in range(5, 25)] \
  + ["phq_%d" % i for i in range(36, 45)] \
  + ["gad_%d" % i for i in range(45, 52)]


def _read_base(path, cohort_label):
    """读取官方基底并数值化，返回分析用 DataFrame（口径与 R read_base 对齐）。"""
    if not os.path.exists(path):
        return None
    d = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    for c in d.columns:
        d[c] = d[c].replace("", np.nan)
    for c in _NUM_COLS:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["cohort"] = cohort_label
    return d


# ---------- 分型（唯一真值，移植自 00_shared_rules.R add_alpha_typing）----------
# 每个 α 各自按本队列 P33 / P67 分「低 / 中 / 高」
#   x <= P33 → 低 ；x >= P67 → 高 ；否则 → 中
# profile = α1cat-α2cat-α3cat-α4cat
# 四者同档 → 归类，否则「其他」
#   低-低-低-低 → 耗竭型 ｜ 中-中-中-中 → 均衡型 ｜ 高-高-高-高 → 高意识型
# self_ok     = phq9_total < 10
# group_type  = 耗竭型&正常→潜在抑郁 / 耗竭型&异常→典型抑郁 / 高意识型&正常→真正健康 / 其余 NA
# potential_depression（潜在痛苦）/ hidden_distress = self_ok & alpha_group=="耗竭型"
# healthy_control = self_ok & alpha_group=="高意识型"
def add_alpha_typing(d):
    for a in ALPHA_AC:
        if a not in d.columns:
            continue
        q = np.nanquantile(d[a].values, [1 / 3, 2 / 3])
        categories = pd.Series(np.where(d[a] <= q[0], "低",
                                np.where(d[a] >= q[1], "高", "中")), index=d.index, dtype=object)
        d[a + "_cat"] = categories.mask(d[a].isna())
    d["profile"] = (d["alpha1_cat"] + "-" + d["alpha2_cat"] + "-" +
                    d["alpha3_cat"] + "-" + d["alpha4_cat"])
    d["alpha_group"] = np.where(d["profile"] == "低-低-低-低", "耗竭型",
                        np.where(d["profile"] == "中-中-中-中", "均衡型",
                         np.where(d["profile"] == "高-高-高-高", "高意识型", "其他")))
    complete_alpha = d[[a + "_cat" for a in ALPHA_AC]].notna().all(axis=1)
    d.loc[~complete_alpha, "alpha_group"] = pd.NA
    d["self_ok"] = d["phq9_total"].lt(10).where(d["phq9_total"].notna()).astype("boolean")
    exhausted = d["alpha_group"].eq("耗竭型").astype("boolean").mask(d["alpha_group"].isna())
    high_awareness = d["alpha_group"].eq("高意识型").astype("boolean").mask(d["alpha_group"].isna())
    so = d["self_ok"]
    d["group_type"] = pd.Series(pd.NA, index=d.index, dtype=object)
    d.loc[(exhausted & so).fillna(False), "group_type"] = "潜在抑郁"
    d.loc[(exhausted & ~so).fillna(False), "group_type"] = "典型抑郁"
    d.loc[(high_awareness & so).fillna(False), "group_type"] = "真正健康"
    d["potential_depression"] = so & exhausted
    d["hidden_distress"] = d["potential_depression"]
    d["healthy_control"] = so & high_awareness
    return d


def add_h_terms(d):
    """公开理论定义 H_trend = α1²；不是后台读数或实测纵向变化。"""
    d["H_trend"] = d["alpha1"] ** 2
    return d


def _finalize(d):
    return add_h_terms(add_alpha_typing(d))


def load_prolific(dedup=True):
    d = _read_base(PROLIFIC_F, "Prolific")
    if d is None:
        return None
    if dedup and "pid" in d.columns:
        # 优先保留 full（与 R 同规则）；基底已每人一行，这里是兜底
        if "data_source" in d.columns:
            d = d.sort_values(["pid", "data_source"],
                             key=lambda s: s.map({"full": 1}).fillna(0),
                             ascending=[True, False])
        d = d.drop_duplicates(subset="pid", keep="first")
    return _finalize(d)


def load_jianshu():
    d = _read_base(JIANSHU_F, "见数")
    if d is None:
        return None
    return _finalize(d)


def load_gold():
    """金标准子样本（仅 Prolific 完整版 n=142），含 MAIA/AAQ/DCS/BADS。"""
    d = _read_base(GOLD_F, "Prolific")
    if d is None:
        return None
    return _finalize(d)


def load_both():
    return {"prolific": load_prolific(), "jianshu": load_jianshu()}


def load_cohort(name):
    if name == "prolific":
        return load_prolific()
    if name == "jianshu":
        return load_jianshu()
    if name == "gold":
        return load_gold()
    raise ValueError("unknown cohort: %s" % name)


# =============================================================
# 统计工具
# =============================================================

def auc_delong(y, x):
    """DeLong 法 AUC 及 95% CI（与 R pROC::ci.auc 对齐）。
    x 已按方向处理：若想表达「x 越高→y 越低」，调用方传入 -x。
    返回 (auc, lo, hi)。
    """
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    m = np.isfinite(y) & np.isfinite(x)
    y, x = y[m], x[m]
    pos = x[y == 1]
    neg = x[y == 0]
    m_, n_ = len(pos), len(neg)
    if m_ == 0 or n_ == 0:
        return (np.nan, np.nan, np.nan)
    auc = ((pos[:, None] > neg[None, :]).sum() +
           0.5 * (pos[:, None] == neg[None, :]).sum()) / (m_ * n_)
    V10 = ((neg[None, :] < pos[:, None]).sum(axis=1) +
           0.5 * (neg[None, :] == pos[:, None]).sum(axis=1)) / n_
    V01 = ((pos[None, :] > neg[:, None]).sum(axis=1) +
           0.5 * (pos[None, :] == neg[:, None]).sum(axis=1)) / m_
    S10 = np.cov(V10) if m_ > 1 else 0.0
    S01 = np.cov(V01) if n_ > 1 else 0.0
    var = S10 / m_ + S01 / n_
    se = np.sqrt(var) if var > 0 else 0.0
    z = 1.959963984540054
    return (auc, max(auc - z * se, 0.0), min(auc + z * se, 1.0))


def auc_oriented(y, x):
    """定向 AUC：返回 max(auc, 1-auc)，使结果恒 >= 0.5（快速判别力指标）。"""
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    m = np.isfinite(y) & np.isfinite(x)
    y, x = y[m], x[m]
    if len(np.unique(y)) < 2:
        return np.nan
    from sklearn.metrics import roc_auc_score
    a = roc_auc_score(y, x)
    return max(a, 1 - a)


def discrimination(y, score):
    """判别力 AUC：score 越高代表越健康（y=1 概率越低）→ 用 -score。
    返回 (auc, lo, hi)，AUC>0.5 表示『score 高→结局少』。"""
    return auc_delong(y, -np.asarray(score, dtype=float))


def cohens_d(x, g):
    x = np.asarray(x, dtype=float)
    g = np.asarray(g)
    a = x[g == 1]
    b = x[g == 0]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    s = np.sqrt(((len(a) - 1) * np.var(a, ddof=1) +
                 (len(b) - 1) * np.var(b, ddof=1)) / (len(a) + len(b) - 2))
    if s == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / s


def pearson_r(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return (np.nan, np.nan)
    r, p = _corr(x, y, method="pearson")
    return (r, p)


def spearman_r(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 3:
        return (np.nan, np.nan)
    from scipy.stats import spearmanr
    r, p = spearmanr(x, y)
    return (r, p)


def _corr(x, y, method="pearson"):
    from scipy.stats import pearsonr
    return pearsonr(x, y)


def cronbach_alpha(items):
    """items: 二维数组 (n_subjects, k_items)，按行计算 Cronbach's α。"""
    items = np.asarray(items, dtype=float)
    items = items[np.all(np.isfinite(items), axis=1)]
    if items.shape[1] < 2:
        return np.nan
    k = items.shape[1]
    item_var = np.var(items, axis=0, ddof=1).sum()
    total_var = np.var(items.sum(axis=1), ddof=1)
    if total_var == 0:
        return np.nan
    return (k / (k - 1)) * (1 - item_var / total_var)


def bootstrap_diff(y1, x1, y2, x2, B=3000, seed=7):
    """固定健康方向的跨队列 AUC 差值及 bootstrap 百分位 CI。
    点估计取原样本 AUC1-AUC2；重抽样不在每次迭代中翻转方向。"""
    rng = np.random.default_rng(seed)
    y1 = np.asarray(y1, dtype=float); x1 = np.asarray(x1, dtype=float)
    y2 = np.asarray(y2, dtype=float); x2 = np.asarray(x2, dtype=float)
    ok1 = np.isfinite(y1) & np.isfinite(x1)
    ok2 = np.isfinite(y2) & np.isfinite(x2)
    y1, x1, y2, x2 = y1[ok1], x1[ok1], y2[ok2], x2[ok2]
    n1, n2 = len(y1), len(y2)
    observed = discrimination(y1, x1)[0] - discrimination(y2, x2)[0]
    ds = []
    for _ in range(B):
        i1 = rng.integers(0, n1, n1)
        i2 = rng.integers(0, n2, n2)
        a1 = discrimination(y1[i1], x1[i1])[0]
        a2 = discrimination(y2[i2], x2[i2])[0]
        if np.isfinite(a1) and np.isfinite(a2):
            ds.append(a1 - a2)
    ds = np.array(ds)
    if len(ds) < 100:
        return (np.nan, np.nan, np.nan, len(ds))
    return (observed, np.quantile(ds, 0.025), np.quantile(ds, 0.975), len(ds))


def logit_or(df, ycol, preds):
    """逻辑回归：返回每预测子的 OR（exp(coef)），基于 statsmodels。
    df 需含 ycol 与 preds 且已 dropna。返回 dict pred->(OR, p)。"""
    import statsmodels.api as sm
    X = df[preds].astype(float).copy()
    X = sm.add_constant(X)
    y = df[ycol].astype(float)
    try:
        # 同一二项最大似然模型使用 GLM/IRLS；避免 Logit/Newton 已到最优解仍误报未收敛。
        res = sm.GLM(y, X, family=sm.families.Binomial()).fit()
        if not bool(res.converged):
            return {p: (np.nan, np.nan) for p in preds}
    except Exception:
        return {p: (np.nan, np.nan) for p in preds}
    out = {}
    for p in preds:
        if p in res.params:
            out[p] = (float(np.exp(res.params[p])), float(res.pvalues[p]))
    return out


def save_tab(df_or_list, name, cohort, module):
    """结果落盘：分析结果/<module>/<cohort>/<name>.csv"""
    d = os.path.join(RESULT_DIR, module, cohort)
    os.makedirs(d, exist_ok=True)
    f = os.path.join(d, name + ".csv")
    if isinstance(df_or_list, (pd.DataFrame, pd.Series)):
        df_or_list.to_csv(f, index=False, encoding="utf-8-sig")
    else:
        pd.DataFrame(df_or_list).to_csv(f, index=False, encoding="utf-8-sig")
    return f


def save_tab_both(df_or_list, name, module):
    """同时写 prolific 与 jianshu 两个子目录（当同一张表同时含两队列时）。"""
    out = []
    for c in COHORTS:
        out.append(save_tab(df_or_list, name, c, module))
    return out


if __name__ == "__main__":
    p = load_prolific()
    j = load_jianshu()
    g = load_gold()
    print("版本:", VERSION_TAG)
    print("Prolific n =", len(p), "| 见数 n =", len(j), "| Gold n =", len(g))
    print("Prolific α 分型分布:")
    print(p["alpha_group"].value_counts(dropna=False).to_dict())
    print("Prolific 潜在痛苦(耗竭型&PHQ<10):", int(p["potential_depression"].sum()))
    print("Prolific 真正健康(高意识型&PHQ<10):", int(p["healthy_control"].sum()))
