#!/usr/bin/env python3
"""Verify example outputs against the bundled reference results."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ASSAYS = ("loMNase", "DNase", "ATAC")
APEX_NUMERIC = (
    "n",
    "apex_x",
    "apex_y_channel_width",
    "pi_enrichment",
    "enrichment_fold",
    "LR",
)


def compare_table(actual_path, expected_path, key, numeric_columns):
    actual = pd.read_csv(actual_path, sep="\t" if actual_path.suffix == ".tsv" else ",")
    expected = pd.read_csv(expected_path, sep="\t" if expected_path.suffix == ".tsv" else ",")

    if key not in actual or key not in expected:
        raise AssertionError(f"missing key column {key!r}")
    actual = actual.sort_values(key).reset_index(drop=True)
    expected = expected.sort_values(key).reset_index(drop=True)
    if actual[key].astype(str).tolist() != expected[key].astype(str).tolist():
        raise AssertionError(f"{key} values differ")

    for column in ("status", "TF_higher", "significant", "gate_passed"):
        if column not in expected:
            continue
        if column not in actual:
            raise AssertionError(f"missing exact-match column {column!r}")
        if actual[column].astype(str).tolist() != expected[column].astype(str).tolist():
            raise AssertionError(f"exact-match column {column!r} differs")

    for column in numeric_columns:
        if column not in expected:
            continue
        if column not in actual:
            raise AssertionError(f"missing numeric column {column!r}")
        observed = pd.to_numeric(actual[column], errors="coerce").to_numpy(float)
        reference = pd.to_numeric(expected[column], errors="coerce").to_numpy(float)
        if not np.allclose(observed, reference, rtol=0.02, atol=0.10, equal_nan=True):
            raise AssertionError(f"numeric column {column!r} differs beyond tolerance")


def main():
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--actual", type=Path, default=here / "results")
    parser.add_argument("--expected", type=Path, default=here / "expected")
    args = parser.parse_args()

    failures = []
    for assay in ASSAYS:
        prefix = f"{assay}_example"
        actual_dir = args.actual / assay
        expected_dir = args.expected / assay
        checks = (
            (f"{prefix}_TF_apex.tsv", "motif", APEX_NUMERIC),
            (f"{prefix}_bias_apex.tsv", "motif", APEX_NUMERIC),
            (
                f"{prefix}_rank_sum_test.csv",
                "feature",
                ("TF_n", "bias_n", "median_TF", "median_bias", "WRS_p", "AUC"),
            ),
            (
                f"{prefix}_cutoffs.csv",
                "feature",
                ("cutoff", "margin", "auc", "loo_acc", "loo_sens", "loo_spec"),
            ),
        )
        for filename, key, numeric in checks:
            actual_path = actual_dir / filename
            expected_path = expected_dir / filename
            try:
                if not actual_path.is_file() or not expected_path.is_file():
                    raise AssertionError("file is missing")
                compare_table(actual_path, expected_path, key, numeric)
                print(f"PASS {assay}: {filename}")
            except Exception as exc:
                failures.append(f"{assay}/{filename}: {exc}")

        for suffix in ("conclusion.txt", "scatter.png", "scatter.pdf"):
            filename = f"{prefix}_{suffix}"
            if not (actual_dir / filename).is_file():
                failures.append(f"{assay}/{filename}: actual file is missing")
            if not (expected_dir / filename).is_file():
                failures.append(f"{assay}/{filename}: expected file is missing")

    if failures:
        print("\nVerification failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print("\nAll three example runs match the bundled expected results.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
