import argparse
import os
import re
import time
import numpy as np
import pandas as pd

from Functions_wrapped import calculate_full_metrics_hierarchical_fast


METRICS_FILE_RE = re.compile(r"^Metrics_N_(?P<N>\d+)_diff_(?P<diff>\d+)\.csv$")


def natural_sort_key(text: str):
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", text)]


def parse_n_diff_from_filename(path: str):
    fname = os.path.basename(path)
    match = METRICS_FILE_RE.match(fname)
    if not match:
        return None, None
    return int(match.group("N")), int(match.group("diff"))


def parse_n_diff_from_dataframe(df: pd.DataFrame):
    if "N" not in df.columns or "diff" not in df.columns or df.empty:
        return None, None

    n_val = pd.to_numeric(df["N"], errors="coerce").dropna()
    diff_val = pd.to_numeric(df["diff"], errors="coerce").dropna()

    if n_val.empty or diff_val.empty:
        return None, None

    return int(n_val.iloc[0]), int(diff_val.iloc[0])


def build_hierarchical_updates(n_val: int, diff_val: int, checks: int):
    hier_metrics = calculate_full_metrics_hierarchical_fast(n_val, diff_val, checks=checks)

    return {
        "N": n_val,
        "diff": diff_val,
        "Mean experiments": float(np.round(hier_metrics[0], 2)),
        "Max samples per pool": int(hier_metrics[1]),
        "N pools": hier_metrics[2],
        "Percentage check": int(hier_metrics[3]),
        "Max experiments per sample": len(hier_metrics[2]) + 1,
        "Mean extra experiments": float(np.round(hier_metrics[4], 2)),
        "Mean steps": float(np.round(hier_metrics[5],2)),
        "Mean positive pools": float(np.round(hier_metrics[6],2)),
        "Mean negative pools": float(np.round(hier_metrics[7],2)),
        "Pooling strategy": "Hierarchical",
    }


def update_metrics_file(path: str, checks: int, inplace: bool = True):
    try:
        df = pd.read_csv(path)
    except Exception as exc:
        print(f"Skipping {path}: could not read CSV ({exc})")
        return False

    if "Pooling strategy" not in df.columns:
        print(f"Skipping {path}: missing 'Pooling strategy' column")
        return False

    n_val, diff_val = parse_n_diff_from_filename(path)
    if n_val is None or diff_val is None:
        n_val, diff_val = parse_n_diff_from_dataframe(df)

    if n_val is None or diff_val is None:
        print(f"Skipping {path}: could not infer N/diff")
        return False

    mask = df["Pooling strategy"].astype(str) == "Hierarchical"
    if not mask.any():
        print(f"Skipping {path}: no Hierarchical row found")
        return False

    try:
        updates = build_hierarchical_updates(n_val, diff_val, checks)
    except Exception as exc:
        print(f"Skipping {path}: failed to recompute Hierarchical metrics ({exc})")
        return False

    hier_idx = df.index[mask][0]

    for col in df.columns:
        if col in updates:
            df.at[hier_idx, col] = updates[col]

    if inplace:
        df.to_csv(path, index=False)
        print(f"Updated Hierarchical row in {path}")

    return True


def main():
    parser = argparse.ArgumentParser(
        description="Recompute and update only Hierarchical rows in Metrics CSV files recursively."
    )
    parser.add_argument(
        "--directory",
        type=str,
        required=True,
        help="Root folder to recursively scan for Metrics_N_*_diff_*.csv files",
    )
    parser.add_argument(
        "--checks",
        type=int,
        default=65,
        help="Checks parameter passed to calculate_full_metrics_hierarchical_fast",
    )
    args = parser.parse_args()

    start_time = time.time()
    scanned = 0
    updated = 0

    for root, dirs, files in os.walk(args.directory, topdown=True):
        dirs[:] = sorted(dirs, key=natural_sort_key)
        for fname in sorted(files, key=natural_sort_key):
            if not METRICS_FILE_RE.match(fname):
                continue

            scanned += 1
            path = os.path.join(root, fname)
            if update_metrics_file(path, checks=args.checks, inplace=True):
                updated += 1

    elapsed = np.round(time.time() - start_time, 2)
    print("\n")
    print("**********************************************************************************************************")
    print(f"Scanned {scanned} metrics files")
    print(f"Updated {updated} Hierarchical rows")
    print(f"Elapsed time: {elapsed} seconds")
    print("**********************************************************************************************************")


if __name__ == "__main__":
    main()
