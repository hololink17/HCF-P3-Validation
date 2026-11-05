# -*- coding: utf-8 -*-
"""Export public L2 analysis bases using an explicit column allowlist.

Only the nine backend-produced readings, external outcomes, demographics, and
version identifiers are eligible for export. Item responses and proprietary
intermediate fields are never selected. Controlled source data are read-only.

Example:
  python export_public_l2.py --source-root /path/to/controlled/project \
      --out-dir /path/to/03_GitHub代码包/code/data/derived
"""
import argparse
import os
import tempfile

import pandas as pd

DIRECT_READINGS = [
    "alpha1", "alpha2", "alpha3", "alpha4",
    "carrier_a", "carrier_b", "carrier_c", "carrier_d", "H_score",
]
PUBLIC_COLUMNS = [
    "participant_token", "cohort", "data_source",
    "questionnaire_version", "mapping_version", "algo_param_version",
    "age_group", "gender", "education",
] + DIRECT_READINGS + [
    "phq9_total", "phq9_normal", "gad7_total", "gad7_item_count",
    "depression", "ibs_positive", "ibs_symptom_count",
    "chronic_pain", "heart_disease", "hypertension", "diabetes",
    "hyperlipidemia", "bmi", "obese", "sleep_hour", "smoke",
    "alcohol_freq", "alcohol_binary", "med_hypertension", "med_diabetes",
    "med_lipid", "med_psych", "sleep_disorder_positive",
    "sleep_symptom_count", "maia_noticing", "maia_not_distracting",
    "maia_not_worrying", "maia_three_dim", "aaq_total", "dcs_total",
    "bads_total",
]
EXPECTED = {"prolific": 462, "jianshu": 504, "gold": 142}
PHASE2_OUTCOMES = ["phq9_total", "ibs_positive", "ibs_symptom_count"]
PHASE2_LONG_COLUMNS = [
    "participant_token", "batch", "order_group", "questionnaire",
    "questionnaire_version", "mapping_version", "algo_param_version", "status",
] + DIRECT_READINGS + PHASE2_OUTCOMES
PHASE2_WIDE_COLUMNS = ["participant_token", "batch", "order_group"] + [
    "%s_%s" % (col, wave)
    for col in DIRECT_READINGS
    for wave in ("new", "old")
]


def _validate_tokens(values, pattern, label, expected_unique=None):
    tokens = pd.Series(list(values), dtype="string")
    if tokens.isna().any() or tokens.eq("").any() or not tokens.str.fullmatch(pattern).all():
        raise AssertionError("%s contains a blank or non-anonymous token" % label)
    if expected_unique is not None and tokens.nunique() != expected_unique:
        raise AssertionError("%s has an invalid unique-token count" % label)
    return tokens


def _atomic_write_outputs(outputs, out_dir):
    """Stage and read back every CSV, then replace destinations with rollback."""
    os.makedirs(out_dir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".p3_l2_stage_", dir=out_dir) as stage_dir:
        staged = {}
        for cohort, frame in outputs.items():
            filename = "p3_base_%s.csv" % cohort
            stage_path = os.path.join(stage_dir, filename)
            frame.to_csv(stage_path, index=False, encoding="utf-8-sig")
            check = pd.read_csv(stage_path, dtype=str, keep_default_na=False,
                                encoding="utf-8-sig")
            expected = frame.astype(str).where(frame.notna(), "").reset_index(drop=True)
            pd.testing.assert_frame_equal(check, expected, check_dtype=False)
            staged[filename] = stage_path

        committed = []
        try:
            for filename, stage_path in staged.items():
                destination = os.path.join(out_dir, filename)
                backup = None
                if os.path.exists(destination):
                    backup = os.path.join(stage_dir, filename + ".before")
                    os.replace(destination, backup)
                committed.append((destination, backup))
                os.replace(stage_path, destination)
        except Exception:
            for destination, backup in reversed(committed):
                if backup is not None:
                    os.replace(backup, destination)
                elif os.path.exists(destination):
                    os.remove(destination)
            raise
        for cohort in outputs:
            print("%s: n=%d, columns=%d -> %s" %
                  (cohort, len(outputs[cohort]), len(outputs[cohort].columns),
                   os.path.join(out_dir, "p3_base_%s.csv" % cohort)))


def _rename_controlled_reading(d, controlled, public):
    out = d.copy()
    if controlled in out.columns:
        if public in out.columns:
            if not out[controlled].equals(out[public]):
                raise AssertionError("controlled/public reading values differ: %s" % public)
            out = out.drop(columns=[controlled])
        else:
            out = out.rename(columns={controlled: public})
        if not out[public].equals(d[controlled]):
            raise AssertionError("reading values changed during rename: %s" % public)
    if public not in out.columns:
        raise AssertionError("phase2 base lacks direct reading: %s" % public)
    return out


def phase2_token_map(long):
    """Build an in-memory crosswalk only when controlled source IDs are present."""
    if "external_profile_id" not in long.columns:
        return None
    ids = long["external_profile_id"].astype(str)
    if ids.eq("").any():
        raise AssertionError("blank controlled phase2 participant ID")
    if "participant_token" in long.columns:
        tokens = long["participant_token"].astype(str)
        if tokens.eq("").any():
            raise AssertionError("blank public phase2 participant token")
        pairs = pd.DataFrame({"id": ids, "token": tokens}).drop_duplicates()
        if pairs["id"].duplicated().any() or pairs["token"].duplicated().any():
            raise AssertionError("phase2 ID/token crosswalk is not one-to-one")
        return dict(zip(pairs["id"], pairs["token"]))
    unique_ids = pd.unique(ids)
    return {pid: "W%04d" % i for i, pid in enumerate(unique_ids, start=1)}


def _attach_phase2_token(d, token_by_id):
    out = d.copy()
    if "participant_token" in out.columns:
        _validate_tokens(out["participant_token"], r"W\d{4}", "phase2 tokens")
        return out
    if token_by_id is None or "external_profile_id" not in out.columns:
        raise AssertionError("controlled phase2 input requires an ID/token crosswalk")
    tokens = out["external_profile_id"].astype(str).map(token_by_id)
    if tokens.isna().any():
        raise AssertionError("phase2 participant ID missing from token crosswalk")
    out.insert(0, "participant_token", tokens)
    return out


def _validate_phase2(long, wide):
    if len(long) != 572 or len(wide) != 166:
        raise AssertionError("phase2 dimensions must be long=572 and wide=166")
    batch = long["batch"].astype(str)
    questionnaire = long["questionnaire"]
    if ((batch == "1") & (questionnaire == "新20")).sum() != 166:
        raise AssertionError("phase2 batch 1 new questionnaire count != 166")
    if ((batch == "1") & (questionnaire == "旧20")).sum() != 166:
        raise AssertionError("phase2 batch 1 old questionnaire count != 166")
    if (batch == "2").sum() != 240:
        raise AssertionError("phase2 batch 2 count != 240")
    _validate_tokens(long["participant_token"], r"W\d{4}", "phase2 long tokens", 406)
    _validate_tokens(wide["participant_token"], r"W\d{4}", "phase2 wide tokens", 166)
    if not wide["participant_token"].is_unique:
        raise AssertionError("phase2 wide tokens are not unique")
    b1_tokens = set(long.loc[batch == "1", "participant_token"])
    if len(b1_tokens) != 166 or b1_tokens != set(wide["participant_token"]):
        raise AssertionError("phase2 tokens differ between batch-1 long and wide")


def sanitize_phase2_long(d, token_by_id=None):
    out = _rename_controlled_reading(d, "h3_present", "H_score")
    out = _attach_phase2_token(out, token_by_id)
    out = out.reindex(columns=[c for c in PHASE2_LONG_COLUMNS if c in out.columns])
    return out


def sanitize_phase2_wide(d, token_by_id=None):
    out = d
    for wave in ("new", "old"):
        out = _rename_controlled_reading(out, "h3_present_" + wave, "H_score_" + wave)
    out = _attach_phase2_token(out, token_by_id)
    out = out.reindex(columns=[c for c in PHASE2_WIDE_COLUMNS if c in out.columns])
    return out


def _read(source_root, cohort):
    path = os.path.join(source_root, "代码", "data", "derived", "p3_base_%s.csv" % cohort)
    d = pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")
    if len(d) != EXPECTED[cohort]:
        raise AssertionError("%s row count %d != %d" % (cohort, len(d), EXPECTED[cohort]))
    return d


def _sanitize(d, tokens, token_pattern):
    missing = [c for c in DIRECT_READINGS if c not in d.columns]
    if missing:
        raise AssertionError("controlled base lacks direct readings: %s" % missing)
    out = d.copy()
    out.insert(0, "public_token", tokens)
    out = out.drop(columns=["participant_token"], errors="ignore")
    out = out.rename(columns={"public_token": "participant_token"})
    out = out.reindex(columns=[c for c in PUBLIC_COLUMNS if c in out.columns])
    _validate_tokens(out["participant_token"], token_pattern, "public participant tokens", len(out))
    if len(out.columns) > len(PUBLIC_COLUMNS):
        raise AssertionError("public allowlist violated")
    return out


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", default=os.environ.get("P3_ROOT"))
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    if not args.source_root:
        raise SystemExit("Specify --source-root or P3_ROOT for the controlled project")
    source_root = os.path.realpath(os.path.abspath(args.source_root))
    out_dir = os.path.realpath(os.path.abspath(args.out_dir))
    controlled_derived = os.path.realpath(os.path.join(source_root, "代码", "data", "derived"))
    if os.path.commonpath([controlled_derived, out_dir]) == controlled_derived:
        raise SystemExit("--out-dir must not be the controlled derived directory or its child")

    prolific = _read(source_root, "prolific")
    jianshu = _read(source_root, "jianshu")
    gold = _read(source_root, "gold")
    p_tokens = ["P%04d" % i for i in range(1, len(prolific) + 1)]
    j_tokens = ["C%04d" % i for i in range(1, len(jianshu) + 1)]
    token_by_pid = dict(zip(prolific["pid"].astype(str), p_tokens))
    g_tokens = [token_by_pid[x] for x in gold["pid"].astype(str)]

    outputs = {
        "prolific": _sanitize(prolific, p_tokens, r"P\d{4}"),
        "jianshu": _sanitize(jianshu, j_tokens, r"C\d{4}"),
        "gold": _sanitize(gold, g_tokens, r"P\d{4}"),
    }

    derived = os.path.join(source_root, "代码", "data", "derived")
    phase2_long = pd.read_csv(os.path.join(derived, "p3_base_phase2_long.csv"),
                              dtype=str, keep_default_na=False, encoding="utf-8-sig")
    phase2_wide = pd.read_csv(os.path.join(derived, "p3_base_phase2_wide.csv"),
                              dtype=str, keep_default_na=False, encoding="utf-8-sig")
    token_by_id = phase2_token_map(phase2_long)
    public_long = sanitize_phase2_long(phase2_long, token_by_id)
    public_wide = sanitize_phase2_wide(phase2_wide, token_by_id)
    _validate_phase2(public_long, public_wide)
    if list(public_long.columns) != PHASE2_LONG_COLUMNS:
        raise AssertionError("public phase2 long does not match its allowlist")
    if list(public_wide.columns) != PHASE2_WIDE_COLUMNS:
        raise AssertionError("public phase2 wide does not match its allowlist")
    batch = public_long["batch"].astype(str)
    if public_long.loc[batch == "2", PHASE2_OUTCOMES].eq("").any().any():
        raise AssertionError("phase2 batch 2 external outcomes are incomplete")
    if not public_long.loc[batch == "1", PHASE2_OUTCOMES].eq("").all().all():
        raise AssertionError("phase2 batch 1 external outcomes must remain missing")
    for col in DIRECT_READINGS:
        source_col = "h3_present" if col == "H_score" else col
        if not public_long[col].equals(phase2_long[source_col]):
            raise AssertionError("phase2 long reading changed: %s" % col)
        for wave in ("new", "old"):
            if not public_wide[col + "_" + wave].equals(
                    phase2_wide[source_col + "_" + wave]):
                raise AssertionError("phase2 wide reading changed: %s_%s" % (col, wave))
    if "external_profile_id" in phase2_long.columns:
        source_ids = set(phase2_long["external_profile_id"].astype(str))
        public_tokens = set(public_long["participant_token"]) | set(public_wide["participant_token"])
        if source_ids & public_tokens:
            raise AssertionError("public phase2 token leaks a controlled participant ID")
        generated = set(token_by_id.values())
        if not all(pd.Series(list(generated)).str.fullmatch(r"W\d{4}")):
            raise AssertionError("generated phase2 tokens are not sequential W tokens")
    outputs["phase2_long"] = public_long
    outputs["phase2_wide"] = public_wide

    from e2_icc_dual_20261031 import build_public_l2 as build_e2_public_l2
    from e2_icc_dual_20261031 import PUBLIC_COLUMNS as E2_PUBLIC_COLUMNS
    public_e2 = build_e2_public_l2(source_root)
    if list(public_e2.columns) != E2_PUBLIC_COLUMNS or len(public_e2) != 164:
        raise AssertionError("public E2 output does not match its allowlist")
    intervals = pd.to_numeric(public_e2["interval_days"], errors="raise")
    if not intervals.between(14, 28).all():
        raise AssertionError("public E2 interval days violate the 14-28 day window")
    outputs["e2_retest"] = public_e2

    shared = sorted(set(outputs["prolific"].columns).intersection(outputs["gold"].columns))
    if not {"participant_token", "phq9_total", "bads_total"}.issubset(shared):
        raise AssertionError("gold/prolific key or outcome columns missing")
    shared = ["participant_token"] + [c for c in shared if c != "participant_token"]
    gold_tokens = set(outputs["gold"]["participant_token"])
    full_shared = outputs["prolific"].loc[
        outputs["prolific"]["participant_token"].isin(gold_tokens), shared]
    gold_shared = outputs["gold"][shared]
    if len(full_shared) != len(gold_shared):
        raise AssertionError("gold/prolific public token coverage differs")
    full_shared = full_shared.sort_values("participant_token").reset_index(drop=True)
    gold_shared = gold_shared.sort_values("participant_token").reset_index(drop=True)
    pd.testing.assert_frame_equal(gold_shared, full_shared, check_dtype=False)
    print("L2 allowlist and full shared-field gold-link checks passed; all validation preceded writes.")
    _atomic_write_outputs(outputs, out_dir)


if __name__ == "__main__":
    main()
