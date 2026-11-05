# -*- coding: utf-8 -*-
"""E2 retest ICC(3,1)/ICC(2,1) from direct backend readouts.

Public mode (default): read code/data/derived/p3_base_e2_retest.csv. The file
contains 164 anonymous pairs, five direct readouts at each wave, and interval
days. It contains no item responses, controlled identifiers, or backend rules.

Controlled export: export_public_l2.py calls build_public_l2(source_root).
Controlled source files are read-only. The analysis never reconstructs scores.

ICC formulas (two measurements, k=2):
  ICC(3,1) = (MSR - MSE) / (MSR + (k-1)*MSE)
  ICC(2,1) = (MSR - MSE) / (MSR + (k-1)*MSE + k*(MSC-MSE)/n)
"""
import csv
import math
import os
import re
from datetime import datetime
from statistics import median

import pandas as pd

ROOT = os.path.abspath(os.environ.get("P3_ROOT") or
                       os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
PUBLIC_L2 = os.path.join(ROOT, "code", "data", "derived", "p3_base_e2_retest.csv")
OUT_DIR = os.environ.get("P3_E2_OUT") or os.path.join(ROOT, "data_output", "E2_results")
OUT_CSV = os.path.join(OUT_DIR, "E2_icc_dual_20261031.csv")
VARIABLES = ["alpha1", "alpha2", "alpha3", "alpha4", "H"]
PUBLIC_COLUMNS = (["participant_token", "interval_days"] +
                  [v + "_time1" for v in VARIABLES] +
                  [v + "_time2" for v in VARIABLES])
ENCODINGS = ("utf-8-sig", "gbk", "utf-8", "gb18030", "latin-1")


def read_csv(path):
    for enc in ENCODINGS:
        try:
            with open(path, encoding=enc, newline="") as f:
                return list(csv.DictReader(f))
        except UnicodeDecodeError:
            continue
    raise RuntimeError("cannot decode: " + path)


def parse_dt(s):
    s = (s or "").strip()
    for fmt in ("%Y-%m-%d %H:%M:%S", "%Y/%m/%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%d %H:%M", "%Y/%m/%d %H:%M"):
        try:
            return datetime.strptime(s, fmt)
        except ValueError:
            pass
    return None


def bridge(paths):
    """Controlled answer ID -> user ID crosswalk; never written to public L2."""
    mapping = {}
    for path in paths:
        for row in read_csv(path):
            aid, uid = row["作答ID"].strip(), row["用户ID"].strip()
            if aid == "作答ID" and uid == "用户ID":
                continue
            assert aid and uid and aid not in mapping, "作答ID缺失或重复"
            mapping[aid] = uid
    return mapping


def scores(row):
    """Read five finalized backend outputs; do not calculate them."""
    values = {v: float(row["scores_" + v]) for v in VARIABLES[:4]}
    values["H"] = float(row["h3_present"])
    return values


def pearson(x, y):
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(x, y))
    sxx = sum((a - mx) ** 2 for a in x)
    syy = sum((b - my) ** 2 for b in y)
    return sxy / math.sqrt(sxx * syy) if sxx > 0 and syy > 0 else float("nan")


def average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i + 1
        while j < len(order) and values[order[j]] == values[order[i]]:
            j += 1
        rank = (i + j + 1) / 2.0
        for k in range(i, j):
            ranks[order[k]] = rank
        i = j
    return ranks


def icc_dual(x, y):
    """Return ICC(3,1), ICC(2,1), Pearson, Spearman, and mean T2-T1."""
    n = len(x)
    k = 2
    gm = sum(x + y) / (2.0 * n)
    rm = [(x[i] + y[i]) / 2.0 for i in range(n)]
    f1, f2 = sum(x) / n, sum(y) / n
    ss_sub = k * sum((m - gm) ** 2 for m in rm)
    ss_for = n * ((f1 - gm) ** 2 + (f2 - gm) ** 2)
    ss_tot = sum((v - gm) ** 2 for v in x + y)
    ss_res = ss_tot - ss_sub - ss_for
    msr = ss_sub / (n - 1)
    msf = ss_for / (k - 1)
    mse = ss_res / ((n - 1) * (k - 1))
    icc31 = (msr - mse) / (msr + (k - 1) * mse)
    icc21 = (msr - mse) / (msr + (k - 1) * mse + k * (msf - mse) / n)
    r = pearson(x, y)
    rho = pearson(average_ranks(x), average_ranks(y))
    mean_t2_minus_t1 = (sum(y) / n) - (sum(x) / n)
    return icc31, icc21, r, rho, mean_t2_minus_t1


def controlled_pairs(source_root):
    """Build paired direct readouts and interval days from controlled sources."""
    source_root = os.path.abspath(source_root)
    exp1_path = os.path.join(source_root, "见数第一期", "后台计算导出", "见数第一期 第一套题 计算导出.csv")
    exp2_paths = [
        os.path.join(source_root, "见数第2期", "后台计算导出", "phase2_A40_先新后旧_旧题计算结果.csv"),
        os.path.join(source_root, "见数第2期", "后台计算导出", "phase2_A40_先旧后新_旧题计算结果.csv"),
    ]
    raw1_paths = [os.path.join(source_root, "见数第一期", "原始数据", "见数%d.csv" % i)
                  for i in range(1, 8)]
    raw2_paths = [
        os.path.join(source_root, "见数第2期", "原始数据", "见数 第二期 40题 先第二套 66人.csv"),
        os.path.join(source_root, "见数第2期", "原始数据", "见数第二期40题 先第二套再第一套20人.csv"),
        os.path.join(source_root, "见数第2期", "原始数据", "见数第二期40题先第一套后第二套30人.csv"),
        os.path.join(source_root, "见数第2期", "原始数据", "见数第二期40题先第一套后第二套50人.csv"),
    ]
    b1, b2 = bridge(raw1_paths), bridge(raw2_paths)
    exp2 = []
    for path in exp2_paths:
        exp2 += read_csv(path)
    t1, t2, t2_raw = {}, {}, {}
    for row in read_csv(exp1_path):
        t1[b1[row["external_profile_id"].strip()]] = scores(row)
    for row in exp2:
        aid = row["external_profile_id"].strip()
        uid = b2[aid]
        t2[uid] = scores(row)
        t2_raw[uid] = aid
    common = sorted(set(t1) & set(t2))
    raw1_by_aid = {}
    for path in raw1_paths:
        for row in read_csv(path):
            raw1_by_aid[row["作答ID"].strip()] = row
    raw2_by_aid = {}
    for path in raw2_paths:
        for row in read_csv(path):
            raw2_by_aid[row["作答ID"].strip()] = row
    aid1_of = {b1[a]: a for a in b1}
    gaps = []
    for uid in common:
        d1 = parse_dt(raw1_by_aid[aid1_of[uid]].get("结束时间"))
        d2 = parse_dt(raw2_by_aid[t2_raw[uid]].get("结束时间"))
        if d1 is None or d2 is None:
            raise AssertionError("controlled retest timestamp is missing or invalid")
        gaps.append((d2 - d1).total_seconds() / 86400.0)
    if len(common) != 164 or len(t1) != 504 or len(t2) != 166:
        raise AssertionError("unexpected controlled E2 dimensions")
    return common, t1, t2, gaps


def build_public_l2(source_root):
    """Export anonymous paired direct readouts; controlled sources stay read-only."""
    common, t1, t2, gaps = controlled_pairs(source_root)
    rows = []
    for i, uid in enumerate(common, start=1):
        row = {"participant_token": "R%04d" % i, "interval_days": gaps[i - 1]}
        for wave, source in (("time1", t1), ("time2", t2)):
            for variable in VARIABLES:
                row[variable + "_" + wave] = source[uid][variable]
        rows.append(row)
    out = pd.DataFrame(rows).reindex(columns=PUBLIC_COLUMNS)
    if out["participant_token"].nunique() != 164:
        raise AssertionError("public E2 tokens are not unique")
    if not out["participant_token"].str.fullmatch(r"R\d{4}").all():
        raise AssertionError("public E2 token format is invalid")
    return out


def load_public_l2(path):
    d = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if list(d.columns) != PUBLIC_COLUMNS or len(d) != 164:
        raise AssertionError("public E2 schema or row count is invalid")
    if not d["participant_token"].str.fullmatch(r"R\d{4}").all() or not d["participant_token"].is_unique:
        raise AssertionError("public E2 participant tokens are invalid")
    if re.search(r"作答ID|用户ID|external_profile_id|submission_id", ",".join(d.columns), re.I):
        raise AssertionError("controlled identifier in public E2 schema")
    common, t1, t2, gaps = [], {}, {}, []
    for row in d.to_dict("records"):
        uid = row["participant_token"]
        common.append(uid)
        gaps.append(float(row["interval_days"]))
        for wave, target in (("time1", t1), ("time2", t2)):
            target[uid] = {v: float(row[v + "_" + wave]) for v in VARIABLES}
    return common, t1, t2, gaps


def main():
    if os.path.exists(PUBLIC_L2):
        common, t1, t2, gaps = load_public_l2(PUBLIC_L2)
        input_source = ("公开L2匿名配对直接读数；配对与重测窗口由受控原始数据核验；"
                        "所有分数均直接读取，不作重构或复算")
        mode = "公开 L2"
    else:
        source_root = os.environ.get("P3_CONTROLLED_ROOT")
        if not source_root:
            raise FileNotFoundError("缺少公开E2基底；受控模式需设置 P3_CONTROLLED_ROOT: " + PUBLIC_L2)
        common, t1, t2, gaps = controlled_pairs(source_root)
        input_source = ("两波后台导出的直出读数；配对键由受控原始数据桥接；"
                        "所有分数均直接读取，不作重构或复算")
        mode = "受控 L3"
    protocol_target_n = 160
    pairable_n = analysis_n = len(common)
    n_note = ("共同可配ID均用于分析，不按方案邀请160人截断；"
              "实际分析n-方案邀请目标n=%+d" % (analysis_n - protocol_target_n))
    sorted_gaps = sorted(gaps)
    median_gap = float(median(sorted_gaps))
    in_window = sum(1 for gap in gaps if 14.0 <= gap <= 28.0)
    if analysis_n != 164 or in_window != analysis_n:
        raise AssertionError("E2 sample size or retest window is invalid")
    print("E2 重测配对（%s）：n=%d；窗口 %d/%d；中位 %.2f 天" %
          (mode, analysis_n, in_window, analysis_n, median_gap))
    rows = []
    for variable in VARIABLES:
        icc31, icc21, r, rho, md = icc_dual(
            [t1[uid][variable] for uid in common], [t2[uid][variable] for uid in common])
        rows.append(dict(proposition="E2_重测", variable=variable, n=analysis_n,
                         protocol_target_n=protocol_target_n, pairable_n=pairable_n,
                         analysis_n=analysis_n, ICC_3_1=round(icc31, 3),
                         ICC_2_1=round(icc21, 3), Pearson_r=round(r, 3),
                         Spearman_rho=round(rho, 3), mean_T2_minus_T1=round(md, 4),
                         n_note=n_note, input_source=input_source))
        print("  %-8s ICC31=%.4f ICC21=%.4f r=%.4f rho=%.4f diff=%+.4f" %
              (variable, icc31, icc21, r, rho, md))
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(OUT_CSV, "w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "proposition", "variable", "n", "protocol_target_n", "pairable_n", "analysis_n",
            "ICC_3_1", "ICC_2_1", "Pearson_r", "Spearman_rho", "mean_T2_minus_T1",
            "n_note", "input_source"])
        writer.writeheader()
        writer.writerows(rows)
    print("E2 结果已写入：" + OUT_CSV)


if __name__ == "__main__":
    main()
