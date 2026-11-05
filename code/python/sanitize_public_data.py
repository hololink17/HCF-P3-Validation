# -*- coding: utf-8 -*-
"""Apply the public L2 allowlist to released analysis bases in place.

This script touches only public derived bases and their released duplicates.
It never reads or modifies controlled raw data. Unknown columns are discarded
by construction; rerunning the script is idempotent.
"""
import argparse
import os

import pandas as pd

from export_public_l2 import (
    DIRECT_READINGS, PUBLIC_COLUMNS, PHASE2_LONG_COLUMNS, PHASE2_WIDE_COLUMNS,
    phase2_token_map, sanitize_phase2_long, sanitize_phase2_wide, _validate_phase2,
    _atomic_write_outputs, _validate_tokens,
)
from e2_icc_dual_20261031 import PUBLIC_COLUMNS as E2_PUBLIC_COLUMNS

PHASE2_LONG = PHASE2_LONG_COLUMNS
PHASE2_WIDE = PHASE2_WIDE_COLUMNS


def _read(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False, encoding="utf-8-sig")


def _allowlisted(d, path, columns, required):
    missing = [c for c in required if c not in d.columns]
    if missing:
        raise AssertionError("%s lacks required direct readings: %s" % (path, missing))
    return d.reindex(columns=[c for c in columns if c in d.columns])


def sanitize(root):
    derived = os.path.join(root, "code", "data", "derived")
    outputs = {}
    for cohort, pattern in (("prolific", r"P\d{4}"), ("jianshu", r"C\d{4}"),
                            ("gold", r"P\d{4}")):
        path = os.path.join(derived, "p3_base_%s.csv" % cohort)
        frame = _allowlisted(_read(path), path, PUBLIC_COLUMNS, DIRECT_READINGS)
        _validate_tokens(frame["participant_token"], pattern, cohort + " tokens", len(frame))
        outputs[cohort] = frame

    long_path = os.path.join(derived, "p3_base_phase2_long.csv")
    long_source = _read(long_path)
    token_by_id = phase2_token_map(long_source)
    long = sanitize_phase2_long(long_source, token_by_id)
    wide_path = os.path.join(derived, "p3_base_phase2_wide.csv")
    wide = sanitize_phase2_wide(_read(wide_path), token_by_id)
    _validate_phase2(long, wide)
    outputs["phase2_long"] = long
    outputs["phase2_wide"] = wide

    e2_path = os.path.join(derived, "p3_base_e2_retest.csv")
    e2 = _allowlisted(_read(e2_path), e2_path, E2_PUBLIC_COLUMNS, E2_PUBLIC_COLUMNS)
    _validate_tokens(e2["participant_token"], r"R\d{4}", "E2 tokens", 164)
    outputs["e2_retest"] = e2
    _atomic_write_outputs(outputs, derived)

    duplicate_dir = os.path.join(root, "data_output", "E3b_results")
    duplicates = {}
    for name, source in (("phase2_long", long), ("phase2_wide", wide)):
        if os.path.exists(os.path.join(duplicate_dir, "p3_base_%s.csv" % name)):
            duplicates[name] = source
    if duplicates:
        _atomic_write_outputs(duplicates, duplicate_dir)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    default_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    parser.add_argument("--root", default=default_root)
    args = parser.parse_args()
    sanitize(os.path.abspath(args.root))
