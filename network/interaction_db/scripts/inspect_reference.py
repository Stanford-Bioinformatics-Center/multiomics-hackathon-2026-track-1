#!/usr/bin/env python3
"""Profile the internal reference parquet file.

Prints schema/dtypes, row count, unique nodes, score quantiles, the share of
non-integer scores, the distribution of fractional parts, reverse-pair and
self-loop counts, the share of rows where protein1 < protein2, the ID length
distribution and whether the rows are sorted by combined_score descending.

Writes the same information to reports/reference_profile.md.
"""
from __future__ import annotations

import argparse
import io
import os
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# Allow running as a plain script.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from string_network.config import load_config  # noqa: E402


def profile(path: str) -> str:
    """Return a markdown report profiling the parquet file at *path*."""
    out = io.StringIO()

    def line(text: str = "") -> None:
        print(text, file=out)

    pf = pq.ParquetFile(path)
    schema = pf.schema_arrow
    df = pd.read_parquet(path)

    line("# Reference file profile")
    line()
    line(f"- **Path:** `{path}`")
    line(f"- **Rows:** {len(df):,}")
    line(f"- **Columns:** {list(df.columns)}")
    line()

    line("## Schema (Arrow)")
    line()
    line("```")
    line(str(schema))
    line("```")
    line()

    line("## Pandas dtypes")
    line()
    line("```")
    line(str(df.dtypes))
    line("```")
    line()

    # Does the parquet store a pandas index?
    has_index = any(
        f.name.startswith("__index_level_") or f.name == "index"
        for f in schema
    )
    line(f"- **Stores a pandas index:** {has_index}")

    p1 = df["protein1"].astype("string")
    p2 = df["protein2"].astype("string")
    score = df["combined_score"].astype("float64")

    nodes = pd.unique(pd.concat([p1, p2], ignore_index=True))
    line(f"- **Unique proteins (nodes):** {len(nodes):,}")
    line()

    line("## Score statistics")
    line()
    q = score.quantile([0, 0.25, 0.5, 0.75, 0.9, 0.99, 1.0])
    line("| quantile | value |")
    line("|---|---|")
    for k, v in q.items():
        line(f"| {k:.2f} | {v:.4f} |")
    line(f"| mean | {score.mean():.4f} |")
    line()

    # Non-integer / fractional analysis.
    frac = (score - np.floor(score)).round(6)
    non_int = frac > 0
    pct_non_int = 100.0 * non_int.mean()
    line(f"- **Non-integer scores:** {non_int.sum():,} ({pct_non_int:.2f}%)")

    # Are all fractions multiples of 0.1?
    tenths = (score * 10).round().astype("int64")
    is_tenth = np.isclose(score, tenths / 10.0, atol=1e-6)
    line(f"- **All scores multiples of 0.1:** {bool(is_tenth.all())}")
    line()

    line("### Distribution of fractional parts (x.1 .. x.9)")
    line()
    frac_digit = (np.round(frac * 10)).astype("int64")
    dist = frac_digit.value_counts().sort_index()
    line("| fractional digit | count | pct |")
    line("|---|---|---|")
    for d, c in dist.items():
        line(f"| .{d} | {c:,} | {100.0 * c / len(df):.2f}% |")
    line()

    # Undirected structure.
    a = np.minimum(p1.to_numpy(dtype=object), p2.to_numpy(dtype=object))
    b = np.maximum(p1.to_numpy(dtype=object), p2.to_numpy(dtype=object))
    canon = pd.Series(list(zip(a, b)))
    dup_canon = canon.duplicated().sum()
    self_loops = int((p1.to_numpy() == p2.to_numpy()).sum())
    dup_rows = int(df.duplicated().sum())
    frac_p1_lt_p2 = float((p1.to_numpy() < p2.to_numpy()).mean())

    # reverse-pair count: pairs where both (A,B) and (B,A) appear as directed rows
    directed = pd.Series(list(zip(p1.to_numpy(object), p2.to_numpy(object))))
    directed_set = set(directed)
    reverse_pairs = sum(1 for (x, y) in directed_set if (y, x) in directed_set and x != y)

    line("## Undirected structure")
    line()
    line(f"- **Self-loops:** {self_loops:,}")
    line(f"- **Fully duplicate rows:** {dup_rows:,}")
    line(f"- **Duplicate undirected canonical keys:** {dup_canon:,}")
    line(f"- **Directed reverse-pairs (A-B and B-A both present):** {reverse_pairs:,}")
    line(f"- **Share of rows with protein1 < protein2:** {100.0 * frac_p1_lt_p2:.2f}%")
    line()

    # ID length distribution.
    line("## ID length distribution")
    line()
    lengths = pd.concat([p1, p2], ignore_index=True).str.len()
    ldist = lengths.value_counts().sort_index()
    line("| length | count |")
    line("|---|---|")
    for L, c in ldist.items():
        line(f"| {L} | {c:,} |")
    has_isoform = pd.concat([p1, p2], ignore_index=True).str.contains(r"-\d+$", regex=True).any()
    line(f"- **Any isoform IDs (e.g. `-2`):** {bool(has_isoform)}")
    line()

    # Sort order.
    is_sorted_desc = bool((score.to_numpy()[:-1] >= score.to_numpy()[1:]).all())
    line("## Ordering")
    line()
    line(f"- **Sorted by combined_score descending:** {is_sorted_desc}")
    line()

    return out.getvalue()


def main() -> None:
    ap = argparse.ArgumentParser(description="Profile the reference parquet file.")
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--path", default=None, help="Override reference path.")
    ap.add_argument("--out", default="reports/reference_profile.md")
    args = ap.parse_args()

    cfg = load_config(args.config)
    path = args.path or cfg["paths"]["reference"]

    report = profile(path)
    print(report)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(report)
    print(f"\nWrote {args.out}")


if __name__ == "__main__":
    main()
