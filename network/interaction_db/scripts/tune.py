#!/usr/bin/env python3
"""Grid-search the settings that matter for matching the reference now that the
score rule is known, and rank each run by match to the reference.

Grid (small and targeted):
    score_mode        : raw  |  physical_else_scaled
    mapping           : idmapping-only (no fallback)
                        fallback, no aliases
                        fallback + UniProt_AC
                        fallback + UniProt_AC + Ensembl_UniProt

The large full/physical links files and each mapping are loaded once and reused.
Writes reports/tuning_results.csv and reports/tuning_summary.md.
"""
from __future__ import annotations

import argparse
import itertools
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from string_network import load as L                                   # noqa: E402
from string_network.build import (BuildOptions, build_network,          # noqa: E402
                                  load_links_full)
from string_network.compare import compare                             # noqa: E402
from string_network.config import load_config                          # noqa: E402

SCORE_MODES = ["raw", "physical_else_scaled"]

# (label, mapping_fallback, alias_sources)
MAPPINGS = [
    ("idmapping-only (no fallback)", False, ()),
    ("fallback, no aliases", True, ()),
    ("fallback +UniProt_AC", True, ("UniProt_AC",)),
    ("fallback +UniProt_AC +Ensembl_UniProt", True, ("UniProt_AC", "Ensembl_UniProt")),
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    args = ap.parse_args()

    cfg = load_config(args.config)
    ref = pd.read_parquet(cfg["paths"]["reference"])
    taxon = int(cfg.get("taxon", 9606))
    threshold = int(cfg.get("threshold", 700))
    scale = float(cfg.get("scale_factor", 0.9))

    print("[tune] loading full + physical links once ...")
    full_links = load_links_full(cfg, links_source="full")
    phys_links = load_links_full(cfg, links_source="physical")
    print(f"[tune] full={len(full_links):,} physical={len(phys_links):,}")

    # Precompute the four mappings.
    mapping_cache = {}
    for label, fb, aliases in MAPPINGS:
        if not fb:
            mapping_cache[label] = L.load_mapping(cfg, "uniprot_idmapping")
        else:
            mapping_cache[label] = L.load_mapping_with_fallback(
                cfg, "uniprot_idmapping", taxon, alias_sources=list(aliases))
        print(f"[tune] mapping '{label}': {len(mapping_cache[label]):,} rows")

    rows = []
    combos = list(itertools.product(SCORE_MODES, MAPPINGS))
    print(f"[tune] evaluating {len(combos)} configurations")
    for i, (score_mode, (label, fb, aliases)) in enumerate(combos, 1):
        opts = BuildOptions(
            threshold=threshold, threshold_when="after", agg="max", decimals=1,
            mapping_source="uniprot_idmapping",
            mapping_fallback=fb, mapping_alias_sources=tuple(aliases),
            links_source="full", score_mode=score_mode, scale_factor=scale,
            taxon=taxon, chunksize=int(cfg.get("chunksize", 2_000_000)),
        )
        t0 = time.time()
        build_df, _ = build_network(
            cfg, opts, links_cache=full_links,
            mapping_cache=mapping_cache[label], physical_cache=phys_links,
            verbose=False)
        mt = compare(build_df, ref)
        dt = time.time() - t0
        rows.append({
            "score_mode": score_mode, "mapping": label,
            "mapping_fallback": fb, "alias_sources": "+".join(aliases) or "none",
            **mt, "seconds": round(dt, 1),
        })
        print(f"[tune] {i}/{len(combos)} score_mode={score_mode:20s} map='{label}' "
              f"| edgeJ={mt['edge_jaccard']:.4f} nodeJ={mt['node_jaccard']:.4f} "
              f"exact={mt['pct_scores_exact']:.1f}% nonint={mt['build_pct_non_int']:.1f}% "
              f"({dt:.1f}s)")

    res = pd.DataFrame(rows)
    reports = cfg["paths"]["reports"]
    os.makedirs(reports, exist_ok=True)
    res.to_csv(os.path.join(reports, "tuning_results.csv"), index=False)

    ranked = res.sort_values(
        ["edge_jaccard", "pct_scores_exact", "node_jaccard"], ascending=False
    ).reset_index(drop=True)
    best = ranked.iloc[0]

    lines = ["# Tuning summary", ""]
    lines.append("Grid over the settings that matter once the score rule is known: "
                 "`score_mode` x mapping fallback / alias sources. Ranked by edge "
                 "Jaccard, then exact-score share, then node Jaccard.")
    lines.append("")
    lines.append("## Best configuration")
    lines.append("")
    lines.append(f"- **score_mode:** `{best['score_mode']}`")
    lines.append(f"- **mapping:** {best['mapping']}")
    lines.append(f"- **mapping_fallback:** `{best['mapping_fallback']}`  "
                 f"**alias_sources:** `{best['alias_sources']}`")
    lines.append("")
    lines.append(f"- Edge Jaccard: **{best['edge_jaccard']:.4f}** "
                 f"(shared {int(best['shared_edges']):,}, "
                 f"only-ref {int(best['edges_only_in_ref']):,}, "
                 f"only-build {int(best['edges_only_in_build']):,})")
    lines.append(f"- Node Jaccard: **{best['node_jaccard']:.4f}** "
                 f"(only-ref {int(best['nodes_only_in_ref']):,}, "
                 f"only-build {int(best['nodes_only_in_build']):,})")
    lines.append(f"- Exact score match on shared edges: {best['pct_scores_exact']:.2f}%")
    lines.append(f"- Non-integer scores: build {best['build_pct_non_int']:.2f}% "
                 f"vs reference {best['ref_pct_non_int']:.2f}%")
    lines.append("")
    lines.append("## All configurations")
    lines.append("")
    cols = ["score_mode", "mapping", "edge_jaccard", "node_jaccard",
            "pct_scores_exact", "build_pct_non_int", "shared_edges",
            "edges_only_in_build", "edges_only_in_ref"]
    lines.append("| " + " | ".join(cols) + " |")
    lines.append("|" + "|".join(["---"] * len(cols)) + "|")
    for _, r in ranked[cols].iterrows():
        vals = []
        for c in cols:
            v = r[c]
            if isinstance(v, float):
                vals.append(f"{v:.4f}" if v < 10 else f"{v:,.2f}")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    lines.append("")

    with open(os.path.join(reports, "tuning_summary.md"), "w") as fh:
        fh.write("\n".join(lines))
    print(f"[tune] wrote {os.path.join(reports, 'tuning_summary.md')}")
    print("BEST:", best["score_mode"], "|", best["mapping"],
          "| edgeJ=%.4f exact=%.2f%%" % (best["edge_jaccard"], best["pct_scores_exact"]))


if __name__ == "__main__":
    main()
