# run_all.py —— P3 工具验证 Python 移植版总驱动
# 仅执行下列公开白名单模块，依次调用 main()，
# 对 Prolific(462) 与 见数第一期(504) 两队列跑公开下游分析。
# 不使用目录自动发现，避免本机保留的受控/历史脚本误入公开回归。
#
# 数据：code/data/derived/ 下六份公开 L2 基底。
# 用法：python run_all.py
import os
import sys
import importlib.util
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

PUBLIC_MODULE_NAMES = [
    "a01_descriptives.py", "a02_alpha_structure.py", "a03_causal_chain.py",
    "a04_predictive.py", "a05_ibs.py", "a06_clinical_utility.py",
    "a07_hidden_distress.py", "a09_carrier.py", "a10_robustness.py",
    "a11_cross_cohort.py", "a12_convergent.py", "a15_screening_roc.py",
    "a16_hierarchical_regression.py", "a18_anchor_vs_algorithm.py",
    "a19_ibs_subgroup.py", "a24_item_equivalence.py",
    "a25_replication_extension.py", "a27_suicidal_ideation.py",
    "a32_design_power.py", "a36_carrier_extended.py",
    "a37_validity_profiles_extended.py", "a38_ibs_clinical_extended.py",
    "a39_anchor_lifestyle_extended.py", "a41_gold_incremental_ipsative.py",
    "a42_semantic_confound.py", "e03b_phase2.py", "e2_icc_dual_20261031.py",
]
MODULES = [os.path.join(HERE, name) for name in PUBLIC_MODULE_NAMES]
missing = [name for name in PUBLIC_MODULE_NAMES if not os.path.isfile(os.path.join(HERE, name))]
if missing:
    raise FileNotFoundError("公开白名单模块缺失：%s" % ", ".join(missing))

PUBLIC_L2_CONTROLLED_ONLY = {
    "a18_anchor_vs_algorithm.py": "需要受控 L3 的 anchor1–anchor4 数值锚点",
    "a27_suicidal_ideation.py": "需要受控 PHQ-9 题项；公开 L2 仅含外部量表总分，受控聚合结果不作为公开可复算结果",
    "a39_anchor_lifestyle_extended.py": "整模块需要受控 L3 的 anchor1–anchor4；公开 L2 只保留未运行状态说明",
}

PUBLIC_L2_PARTIAL = {
    "a02_alpha_structure.py": "公开 L2 不含题目级作答；题目级 PCA 未运行，其余结构分析已完成",
    "a07_hidden_distress.py": "公开 L2 不含 PHQ/GAD 逐题输入；疲乏条目及双正常部分题项级分支未运行，其余九读数与结局分析已完成",
    "a10_robustness.py": "公开 L2 不含注意力字段；注意力敏感性未运行，分层稳健性已完成",
    "a12_convergent.py": "公开 L2 不含时长/注意力字段；可接受性未运行，Gold 收敛分析已完成",
    "a24_item_equivalence.py": "公开 L2 不含题目级输入；题目等价性未运行，仅完成症状分层描述",
    "a25_replication_extension.py": "公开 L2 不含锚点题；锚定信度未运行，其余扩展分析已完成",
    "a36_carrier_extended.py": "稀疏/未收敛等不可估计模型列入错误表；模块为部分完成而非全部模型成功",
}


def public_l2_mode():
    """以公开基底的实际字段识别 L2，避免把受控模块误报为公开可运行。"""
    import pandas as pd
    root = os.path.abspath(os.environ.get("P3_ROOT", os.path.join(HERE, os.pardir, os.pardir)))
    f = os.path.join(root, "code", "data", "derived", "p3_base_prolific.csv")
    if not os.path.exists(f):
        raise FileNotFoundError("缺少分析基底：%s" % f)
    cols = pd.read_csv(f, nrows=0, encoding="utf-8-sig").columns
    controlled_markers = {"pid", "hcf_5", "phq_36", "anchor1"}
    return not any(c in cols for c in controlled_markers)


def run_one(path):
    name = os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # e03b_phase2 is an output-level script whose top-level block is the entrypoint.
    # All aNN modules expose main(); do not execute E3b twice or report a false failure.
    if name == "e03b_phase2":
        return
    mod.main()


def main():
    print("=" * 70)
    print("P3 工具验证 · Python 移植版 · 全模块运行")
    print("版本口径：问卷与映射 科研-第一套题 ｜ 算法参数 算法v1.1 ｜ 单一固定管线")
    print("模块数：%d" % len(MODULES))
    print("=" * 70)
    is_public_l2 = public_l2_mode()
    print("数据模式：%s" % ("公开 L2" if is_public_l2 else "受控 L3 完整基底"))
    ok, failed, skipped, partial = [], [], [], []
    reasons = {}
    for path in MODULES:
        name = os.path.basename(path)
        if is_public_l2 and name in PUBLIC_L2_CONTROLLED_ONLY:
            skipped.append(name)
            reasons[name] = PUBLIC_L2_CONTROLLED_ONLY[name]
            print("\n>>> 跳过 %s：%s" % (name, reasons[name]))
            continue
        print("\n>>> 运行 %s" % name)
        try:
            run_one(path)
            if is_public_l2 and name in PUBLIC_L2_PARTIAL:
                partial.append(name)
                reasons[name] = PUBLIC_L2_PARTIAL[name]
                print(">>> 部分完成 %s：%s" % (name, reasons[name]))
            else:
                ok.append(name)
        except Exception as e:
            failed.append(name)
            reasons[name] = str(e)
            print("!!! %s 失败: %s" % (name, e))
            traceback.print_exc()
    print("\n" + "=" * 70)
    print("完整成功 %d；部分完成 %d；跳过 %d；失败 %d；总计 %d" %
          (len(ok), len(partial), len(skipped), len(failed), len(MODULES)))
    public_root = os.path.abspath(os.environ.get(
        "P3_ROOT", os.path.join(HERE, os.pardir, os.pardir)))
    result_dir = os.environ.get("P3_RESULT_DIR")
    manifest = (os.path.join(os.path.abspath(result_dir), os.pardir, "run_all_manifest.csv")
                if result_dir else os.path.join(public_root, "data_output", "run_all_manifest.csv"))
    os.makedirs(os.path.dirname(manifest), exist_ok=True)
    import pandas as pd
    rows = []
    for p in MODULES:
        name = os.path.basename(p)
        status = ("success" if name in ok else
                  "partial_public_l2" if name in partial else
                  "skipped_controlled_l3_required" if name in skipped else "failed")
        rows.append({"module": name, "status": status, "reason": reasons.get(name, "")})
    pd.DataFrame(rows).to_csv(manifest, index=False, encoding="utf-8-sig")
    if failed:
        print("失败：%s" % ", ".join(failed))
        raise SystemExit(1)
    if partial:
        print("公开 L2 部分完成：%s" % ", ".join(partial))
        print("部分完成不等同于该模块声明的全部子分析均已执行；详见 reason。")
    if skipped:
        print("公开 L2 受控输入跳过：%s" % ", ".join(skipped))
        print("受控 L3 专属模块未执行，不计为成功。")
    if not partial and not skipped:
        print("✓ 全部模块完整成功")


if __name__ == "__main__":
    main()
