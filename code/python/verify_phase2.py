# verify_phase2.py —— 公开 L2 第二期基底完整性与公开边界核验
# 只验证后台直出读数与公开字段，不读取原始作答，也不重构任何内部分数。
import os
from pathlib import Path

import pandas as pd

from sanitize_public_data import PHASE2_LONG, PHASE2_WIDE

ROOT = Path(os.environ.get("P3_ROOT", Path(__file__).resolve().parents[2])).resolve()
DERIVED = ROOT / "code" / "data" / "derived"
SCORES = ["alpha1", "alpha2", "alpha3", "alpha4",
          "carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score"]

long = pd.read_csv(DERIVED / "p3_base_phase2_long.csv", dtype=str,
                   keep_default_na=False, encoding="utf-8-sig")
wide = pd.read_csv(DERIVED / "p3_base_phase2_wide.csv", dtype=str,
                   keep_default_na=False, encoding="utf-8-sig")
assert len(long) == 572 and len(wide) == 166
assert list(long.columns) == PHASE2_LONG, "第二期长表必须精确匹配公开schema及列顺序"
assert list(wide.columns) == PHASE2_WIDE, "第二期宽表必须精确匹配公开schema及列顺序"
assert long["participant_token"].str.fullmatch(r"W\d{4}").all()
assert wide["participant_token"].str.fullmatch(r"W\d{4}").all()
assert long["participant_token"].nunique() == 406
assert wide["participant_token"].nunique() == 166
assert wide["participant_token"].is_unique
assert not long.duplicated(["participant_token", "questionnaire"]).any(), "同一参与者同一问卷重复"
assert set(long["batch"]) == {"1", "2"}
assert set(long["status"]) == {"done"}

expected_groups = {
    ("1", "先新后旧", "新20"): 86, ("1", "先新后旧", "旧20"): 86,
    ("1", "先旧后新", "新20"): 80, ("1", "先旧后新", "旧20"): 80,
    ("2", "单批", "新20"): 240,
}
actual_groups = long.groupby(["batch", "order_group", "questionnaire"]).size().to_dict()
assert actual_groups == expected_groups, "第二期批次/顺序/问卷组合计数不符"
b1 = long.loc[long["batch"] == "1"].copy()
b2 = long.loc[long["batch"] == "2"].copy()
assert set(b1["participant_token"]) == set(wide["participant_token"])
assert b1.groupby("participant_token")["questionnaire"].agg(set).eq({"新20", "旧20"}).all()
assert b2.groupby("participant_token")["questionnaire"].agg(set).eq({"新20"}).all()
wide_meta = wide[["participant_token", "batch", "order_group"]].sort_values("participant_token").reset_index(drop=True)
long_meta = b1[["participant_token", "batch", "order_group"]].drop_duplicates().sort_values(
    "participant_token").reset_index(drop=True)
pd.testing.assert_frame_equal(wide_meta, long_meta)
for frame in (long, wide):
    assert not any("external_profile_id" in c or "submission_id" in c or
                   "submitted_at" in c for c in frame.columns), "公开包不得包含平台ID或时间戳"

for c in SCORES:
    long[c] = pd.to_numeric(long[c], errors="raise")
    assert long[c].notna().all()
    assert long[c].between(0, 1.0001).all(), "%s 超出[0,1]" % c
for wave in ("new", "old"):
    for c in SCORES:
        col = c + "_" + wave
        wide[col] = pd.to_numeric(wide[col], errors="raise")
        assert wide[col].notna().all()
        assert wide[col].between(0, 1.0001).all(), "%s 超出[0,1]" % col
for c in SCORES:
    b1 = long.loc[long["batch"] == "1"]
    paired = b1.pivot(index="participant_token", columns="questionnaire", values=c).rename(
        columns={"新20": c + "_new", "旧20": c + "_old"}).sort_index()
    paired.columns.name = None
    wide_pair = wide.set_index("participant_token")[[c + "_new", c + "_old"]].sort_index()
    pd.testing.assert_frame_equal(paired, wide_pair, check_dtype=False)

assert set(long["algo_param_version"]) == {"算法v1.1"}
for questionnaire, version in (("新20", "科研-第二套题"), ("旧20", "科研-第一套题")):
    rows = long.loc[long["questionnaire"] == questionnaire]
    assert set(rows["questionnaire_version"]) == {version}
    assert set(rows["mapping_version"]) == {version}
assert b2[["phq9_total", "ibs_positive", "ibs_symptom_count"]].ne("").all().all()
assert b1[["phq9_total", "ibs_positive", "ibs_symptom_count"]].eq("").all().all()
phq = pd.to_numeric(b2["phq9_total"], errors="raise")
ibs = pd.to_numeric(b2["ibs_positive"], errors="raise")
ibs_count = pd.to_numeric(b2["ibs_symptom_count"], errors="raise")
assert phq.between(0, 27).all() and phq.mod(1).eq(0).all()
assert ibs.isin([0, 1]).all()
assert ibs_count.ge(0).all() and ibs_count.mod(1).eq(0).all()

print("公开第二期基底核验通过：long=572行/406人，wide=166对。")
print("九个后台直出读数均为数值且范围有效；未执行任何内部计分重构。")
print("\n各组去重人数：")
print(long.groupby(["batch", "order_group", "questionnaire"])["participant_token"].nunique())
print("\nalpha均值/SD：")
print(long.groupby(["batch", "questionnaire"])[SCORES[:4]].agg(["mean", "std"]).round(4))
