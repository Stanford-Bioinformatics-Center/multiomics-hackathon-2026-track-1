"""Compare a built network against the reference and write a report."""
from __future__ import annotations

import os
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd


def _canon_keys(df: pd.DataFrame) -> pd.Series:
    """Undirected canonical 'a\\tb' key per row (a <= b)."""
    p1 = df["protein1"].to_numpy(dtype=object)
    p2 = df["protein2"].to_numpy(dtype=object)
    a = np.where(p1 <= p2, p1, p2)
    b = np.where(p1 <= p2, p2, p1)
    return pd.Series([f"{x}\t{y}" for x, y in zip(a, b)])


def _nodes(df: pd.DataFrame) -> set:
    return set(pd.concat([df["protein1"], df["protein2"]], ignore_index=True).tolist())


def compare(build_df: pd.DataFrame, ref_df: pd.DataFrame) -> Dict[str, float]:
    """Return a dict of comparison metrics between build and reference."""
    bkey = _canon_keys(build_df)
    rkey = _canon_keys(ref_df)

    bmap = dict(zip(bkey, build_df["combined_score"].astype(float).to_numpy()))
    rmap = dict(zip(rkey, ref_df["combined_score"].astype(float).to_numpy()))

    bset, rset = set(bmap), set(rmap)
    inter = bset & rset
    union = bset | rset

    edge_jaccard = len(inter) / len(union) if union else 0.0

    bnodes, rnodes = _nodes(build_df), _nodes(ref_df)
    ninter = bnodes & rnodes
    nunion = bnodes | rnodes
    node_jaccard = len(ninter) / len(nunion) if nunion else 0.0

    # Score agreement on shared edges.
    if inter:
        bs = np.array([bmap[k] for k in inter])
        rs = np.array([rmap[k] for k in inter])
        diff = np.abs(bs - rs)
        pct_exact = 100.0 * np.mean(diff < 0.05)
        mad = float(np.mean(diff))
        corr = float(np.corrcoef(bs, rs)[0, 1]) if len(inter) > 1 else float("nan")
    else:
        pct_exact = mad = corr = float("nan")

    def pct_non_int(scores: np.ndarray) -> float:
        frac = scores - np.floor(scores)
        return 100.0 * np.mean(frac > 1e-6)

    return {
        "ref_edges": len(rset),
        "build_edges": len(bset),
        "shared_edges": len(inter),
        "edges_only_in_ref": len(rset - bset),
        "edges_only_in_build": len(bset - rset),
        "edge_jaccard": edge_jaccard,
        "ref_nodes": len(rnodes),
        "build_nodes": len(bnodes),
        "shared_nodes": len(ninter),
        "nodes_only_in_ref": len(rnodes - bnodes),
        "nodes_only_in_build": len(bnodes - rnodes),
        "node_jaccard": node_jaccard,
        "pct_scores_exact": pct_exact,
        "mean_abs_diff": mad,
        "score_corr": corr,
        "build_pct_non_int": pct_non_int(build_df["combined_score"].astype(float).to_numpy()),
        "ref_pct_non_int": pct_non_int(ref_df["combined_score"].astype(float).to_numpy()),
    }


def format_report(metrics: Dict[str, float], title: str = "Comparison vs reference") -> str:
    def f(x):
        return f"{x:,.4f}" if isinstance(x, float) else f"{x:,}"

    lines = [f"# {title}", ""]
    lines.append("## Edges")
    lines.append("")
    lines.append(f"- Reference edges: {metrics['ref_edges']:,}")
    lines.append(f"- Build edges: {metrics['build_edges']:,}")
    lines.append(f"- Shared edges: {metrics['shared_edges']:,}")
    lines.append(f"- Only in reference: {metrics['edges_only_in_ref']:,}")
    lines.append(f"- Only in build: {metrics['edges_only_in_build']:,}")
    lines.append(f"- **Edge Jaccard: {metrics['edge_jaccard']:.4f}**")
    lines.append("")
    lines.append("## Nodes")
    lines.append("")
    lines.append(f"- Reference nodes: {metrics['ref_nodes']:,}")
    lines.append(f"- Build nodes: {metrics['build_nodes']:,}")
    lines.append(f"- Shared nodes: {metrics['shared_nodes']:,}")
    lines.append(f"- Only in reference: {metrics['nodes_only_in_ref']:,}")
    lines.append(f"- Only in build: {metrics['nodes_only_in_build']:,}")
    lines.append(f"- **Node Jaccard: {metrics['node_jaccard']:.4f}**")
    lines.append("")
    lines.append("## Scores on shared edges")
    lines.append("")
    lines.append(f"- % exactly equal (|diff| < 0.05): {metrics['pct_scores_exact']:.2f}%")
    lines.append(f"- Mean absolute difference: {metrics['mean_abs_diff']:.4f}")
    lines.append(f"- Correlation: {metrics['score_corr']:.4f}")
    lines.append("")
    lines.append("## Non-integer scores")
    lines.append("")
    lines.append(f"- Build: {metrics['build_pct_non_int']:.2f}%")
    lines.append(f"- Reference: {metrics['ref_pct_non_int']:.2f}%")
    lines.append("")
    return "\n".join(lines)


def compare_files(build_path: str, ref_path: str,
                  out_path: Optional[str] = None,
                  title: str = "Comparison vs reference") -> Tuple[Dict[str, float], str]:
    build_df = pd.read_parquet(build_path)
    ref_df = pd.read_parquet(ref_path)
    metrics = compare(build_df, ref_df)
    report = format_report(metrics, title=title)
    if out_path:
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, "w") as fh:
            fh.write(report)
        print(f"[compare] wrote {out_path}")
    return metrics, report
