# verify_scoring.py —— 计分口径核对（针对「计分错误/反向计分」专项核查）
# 核验：金标准公式中反向计分配置、基底内 PHQ-9/GAD-7 总分、MAIA 三维均值、
#      MAIA/AAQ/DCS/BADS 官方计分一致性。输出 分析结果/计分核对/ 报告。
import os
import json
import numpy as np
import pandas as pd
import common as C

# 受控项目使用中文“代码”目录；显式切换数据位置，但不改变公开包默认路径。
CONTROLLED_ROOT = os.environ.get("P3_CONTROLLED_ROOT")
if CONTROLLED_ROOT:
    controlled_derived = os.path.join(os.path.abspath(CONTROLLED_ROOT), "代码", "data", "derived")
    C.DERIVED = controlled_derived
    C.PROLIFIC_F = os.path.join(controlled_derived, "p3_base_prolific.csv")
    C.JIANSHU_F = os.path.join(controlled_derived, "p3_base_jianshu.csv")
    C.GOLD_F = os.path.join(controlled_derived, "p3_base_gold.csv")
    C.PROLIFIC_RAW_CLEAN = os.path.join(
        os.path.abspath(CONTROLLED_ROOT), "Prolific", "原始数据", "prolific_clean_人口学PHQ9疾病.csv")

GOLD_JSON = os.path.join(C.CODE_DIR, "config", "gold_formula.json")
MAIA_JSON = os.path.join(C.CODE_DIR, "config", "maia_formula_v3.json")
OUT = os.path.join(C.RESULT_DIR, "计分核对")
os.makedirs(OUT, exist_ok=True)

# 独立核验直接写死官方原始题项规则，不调用 build_analysis_base.score_gold，
# 也不以 gold_formula.json 作为计分正确性的唯一依据。
OFFICIAL_RAW = {
    "aaq_total": (list(range(66, 73)), set(), 1, 7),
    "dcs_total": (list(range(73, 89)),
                  {75, 76, 77, 78, 80, 81, 84, 85, 86, 87, 88}, 1, 5),
    "bads_total": (list(range(89, 114)),
                   {89, 94, 96, 97, 101, 102, 104, 105, 107, 108, 109, 110, 112, 113}, 0, 6),
}


def write_reports(rows, public_l2):
    """保存实际执行/跳过状态，不把缺少题项写成计分正确。"""
    rep = pd.DataFrame(rows, columns=["类别", "项目", "核对结果"])
    rep.to_csv(os.path.join(OUT, "计分口径核对.csv"), index=False, encoding="utf-8-sig")
    if public_l2:
        verdict_rows = [
            ("运行边界", "公开 L2：只核验公开总分/衍生字段；原始题项重算均标为受控跳过"),
            ("PHQ-9/GAD-7", "仅核验公开总分范围；未核验逐题0–3或总分=Σ题项"),
            ("MAIA-2/AAQ-II/DCS/BADS", "公开原始题项缺失，未执行逐人官方规则重算；配置列出不等于独立验证"),
            ("九个后台直出读数", "只读取直出字段，不重构或验证内部计分"),
            ("疾病/IBS", "仅核验公开二值字段取值；未由公开题项独立核验 Rome IV 结构"),
            ("潜在痛苦(A4)", "本脚本不独立核验分组定义；其实现位于 common.py，不能把本报告当作该定义的验证"),
        ]
    else:
        verdict_rows = [
            ("运行边界", "受控 L3：已执行可用的原始题项重算；仍未验证专有后台内部计算"),
            ("PHQ-9/GAD-7", "逐题0–3且总分=Σ题项；以实际核验表中的n与最大偏差为准"),
            ("MAIA-2/AAQ-II/DCS/BADS", "按当前官方规则配置完成受控逐人重算；以实际核验表为准"),
            ("九个后台直出读数", "只读取直出字段，不重构或验证内部计分"),
            ("疾病/IBS", "本脚本仅核查二值取值；Rome IV 结构须由清洗链另行核验"),
            ("潜在痛苦(A4)", "本脚本不独立核验分组定义；其实现位于 common.py"),
        ]
    verdict = pd.DataFrame(verdict_rows, columns=["维度", "结论"])
    verdict.to_csv(os.path.join(OUT, "计分正确性结论.csv"), index=False, encoding="utf-8-sig")
    print(rep.to_string(index=False))
    print("\n结论：")
    print(verdict.to_string(index=False))


def main():
    rows = []
    # 1. 金标准公式 JSON：反向计分配置
    gold = json.load(open(GOLD_JSON, encoding="utf-8"))
    maia = json.load(open(MAIA_JSON, encoding="utf-8"))
    for k in ["aaq_total", "dcs_total", "bads_total"]:
        rows.append(("金标准公式", k, "正向题数=%d 反向题数=%d max_abs_err=%s"
                     % (gold[k]["n_pos"], gold[k]["n_rev"], gold[k]["max_abs_err"])))
    for k, v in maia.items():
        if "items" in v:
            rev = sum(1 for x in v["items"].values() if x == -1)
            pos = sum(1 for x in v["items"].values() if x == 1)
            rows.append(("MAIA公式", k, "正向=%d 反向=%d R2=%s" % (pos, rev, v.get("R2"))))

    p = C.load_prolific()
    public_l2 = not bool(CONTROLLED_ROOT)
    if not public_l2:
        required = [C.PROLIFIC_F, C.JIANSHU_F, C.GOLD_F, C.PROLIFIC_RAW_CLEAN]
        missing = [f for f in required if not os.path.isfile(f)]
        if missing:
            raise FileNotFoundError("受控计分核验缺少输入：%s" % "; ".join(missing))
    rows.append(("运行模式", "输入层级", "公开 L2" if public_l2 else "受控 L3"))
    # 2/3. 只有受控题项存在时才核验逐题范围与总分求和。
    for cohort, dd in [("prolific", p), ("jianshu", C.load_jianshu())]:
        phq_items = ["phq_%d" % i for i in range(36, 45)]
        gad_items = ["gad_%d" % i for i in range(45, 52)]
        for label, items, total, total_max in [("PHQ-9", phq_items, "phq9_total", 27),
                                               ("GAD-7", gad_items, "gad7_total", 21)]:
            if all(c in dd.columns for c in items + [total]):
                vals = dd[items]
                rng = (vals.min().min(), vals.max().max())
                calc = vals.sum(axis=1, min_count=len(items))
                both = calc.notna() & dd[total].notna()
                if not public_l2 and (not bool(both.all()) or vals.isna().any().any()):
                    raise AssertionError("%s %s 受控逐题核验要求题项与总分完整" % (cohort, label))
                err = (calc[both] - dd.loc[both, total]).abs().max()
                if rng[0] < 0 or rng[1] > 3 or (both.sum() and err > 1e-12):
                    raise AssertionError("%s %s 计分核验失败：range=%s max_err=%s" % (cohort, label, rng, err))
                rows.append(("基底核验", cohort + " " + label + " 总分",
                             "题项0–3；总分=Σ题项；最大偏差=%.8f；完整n=%d" % (err, both.sum())))
            elif total in dd.columns:
                if not public_l2:
                    raise AssertionError("%s %s 受控模式缺少完整逐题列" % (cohort, label))
                values = pd.to_numeric(dd[total], errors="raise")
                assert values.dropna().between(0, total_max).all(), "%s %s 总分越界" % (cohort, label)
                rows.append(("公开字段核验", cohort + " " + label + " 总分",
                             "仅核验总分范围0–%d；公开 L2 无逐题，未核验逐题0–3或总分=Σ题项" % total_max))
            else:
                rows.append(("输入缺失", cohort + " " + label, "缺少题项与总分列"))
    if public_l2:
        for target in ["aaq_total", "dcs_total", "bads_total"]:
            rows.append(("受控核验跳过", target, "公开 L2 无原始题项；未执行官方规则逐人重算"))
        for target, spec in maia.items():
            if "items" in spec:
                rows.append(("受控核验跳过", target, "公开 L2 无原始题项；未执行5-x反向及逐人重算"))
        gold_public = C.load_gold()
        maia_sub = ["maia_noticing", "maia_not_distracting", "maia_not_worrying"]
        if all(c in gold_public.columns for c in maia_sub + ["maia_three_dim"]):
            md_calc = gold_public[maia_sub].dropna().mean(axis=1)
            md_err = (md_calc - gold_public.loc[md_calc.index, "maia_three_dim"]).abs().max()
            assert np.isfinite(md_err) and md_err <= 1e-12
            rows.append(("衍生字段核验", "MAIA 三维均值",
                         "mean(3子量表) 与基底最大偏差=%.8f；不替代原始题项核验" % md_err))
        for c in ["depression", "ibs_positive", "chronic_pain", "hypertension", "diabetes"]:
            if c in p.columns:
                vc = p[c].dropna().unique()
                assert set(vc).issubset({0, 1}), "%s 非二值" % c
                rows.append(("公开二值字段核验", c, "取值集合=%s；不代表已由题项重算" % sorted(vc.tolist())))
        write_reports(rows, True)
        return
    # 4. Gold 四量表：受控 L3 下从原始题项独立重算，与142人Gold基底逐人对照。
    raw = pd.read_csv(C.PROLIFIC_RAW_CLEAN, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    raw = raw.loc[raw["data_source"].eq("full")].copy().set_index("participant_token", drop=False)
    gold_base = C.load_gold().set_index("pid", drop=False)
    assert len(raw) == len(gold_base) == 142 and raw.index.is_unique and gold_base.index.is_unique
    raw = raw.loc[gold_base.index].copy()
    for target in ["aaq_total", "dcs_total", "bads_total"]:
        numbers, reverse, minimum, maximum = OFFICIAL_RAW[target]
        prefix = target.split("_")[0]
        calc = pd.Series(0.0, index=raw.index, dtype=float)
        naive = pd.Series(0.0, index=raw.index, dtype=float)
        for i in numbers:
            item = "%s_%d" % (prefix, i)
            x = pd.to_numeric(raw[item], errors="raise")
            assert x.notna().all() and x.between(minimum, maximum).all(), "%s 原始题项缺失或越界" % item
            calc = calc + (6 - x if i in reverse else x)
            naive = naive + x
        # 用全低/全高合成作答核验反向规则确实生效，而非只复现配置。
        if target == "dcs_total":
            low = sum(6 - minimum if i in reverse else minimum for i in numbers)
            high = sum(6 - maximum if i in reverse else maximum for i in numbers)
            assert (low, high) == (60, 36), (target, low, high)
        elif target == "bads_total":
            low = sum(6 - minimum if i in reverse else minimum for i in numbers)
            high = sum(6 - maximum if i in reverse else maximum for i in numbers)
            assert (low, high) == (84, 66), (target, low, high)
        else:
            assert (7 * minimum, 7 * maximum) == (7, 49)
        if target != "aaq_total":
            assert not np.allclose(calc, naive), "%s 官方计分不应等于未反向裸求和" % target
        err = (calc - gold_base[target]).abs().max()
        if len(calc) != 142 or err > 1e-12:
            raise AssertionError("Gold %s 官方规则独立重算不一致 n=%d err=%s" % (target, len(calc), err))
        rows.append(("Gold官方规则独立重算", target,
                     "142/142；反向题=%d；最大偏差=%.8f" % (len(reverse), err)))
    for target, spec in maia.items():
        if "items" not in spec:
            continue
        weights = np.array(list(spec["items"].values()), dtype=float)
        n_rev = int((weights == -1).sum())
        calc = pd.Series(0.0, index=raw.index)
        for item, weight in spec["items"].items():
            x = pd.to_numeric(raw[item], errors="raise")
            assert x.notna().all() and x.between(0, 5).all(), "%s MAIA原始题项缺失或越界" % item
            if weight == 1:
                calc = calc + x
            elif weight == -1:
                calc = calc + (5 - x)
            else:
                raise AssertionError("MAIA 非单位权重：%s=%s" % (item, weight))
        calc = calc / len(spec["items"])
        err = (calc - gold_base[target]).abs().max()
        if len(calc) != 142 or err > 1e-12:
            raise AssertionError("Gold %s 独立重算不一致 n=%d err=%s" % (target, len(calc), err))
        rows.append(("Gold独立重算", target, "142/142；反向题=%d；最大偏差=%.8f" % (n_rev, err)))
    # 5. MAIA 三维 = mean(觉察,不分散,不担忧) 与基底 maia_three_dim 比较
    maia_sub = ["maia_noticing", "maia_not_distracting", "maia_not_worrying"]
    if all(c in p.columns for c in maia_sub) and "maia_three_dim" in p.columns:
        md_calc = p[maia_sub].dropna().mean(axis=1)
        md_err = (md_calc - p.loc[md_calc.index, "maia_three_dim"]).abs().max()
        assert np.isfinite(md_err) and md_err <= 1e-12, "MAIA 三维均值与基底不一致"
        rows.append(("基底核验", "MAIA 三维均值", "mean(3子量表) 与基底最大偏差=%.8f" % md_err))
    # 6. 疾病/IBS 二值化（"不知道"→0）
    for c in ["depression", "ibs_positive", "chronic_pain", "hypertension", "diabetes"]:
        if c in p.columns:
            vc = p[c].dropna().unique()
            rows.append(("二值化", c, "取值集合=%s" % sorted(vc.tolist())))

    write_reports(rows, False)


if __name__ == "__main__":
    main()
