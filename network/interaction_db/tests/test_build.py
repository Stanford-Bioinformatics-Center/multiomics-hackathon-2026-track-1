"""Tests for the build pipeline using tiny synthetic fixtures.

Fixture design (idmapping.dat, STRING rows only):
    ENSPA -> {P1}          ENSPB -> {P2, P2B}     ENSPC -> {P3}
    ENSPD -> {P4}  (from P4-2, isoform stripped)  ENSPE -> {P5, P1}
    (a non-STRING 'Ensembl' row is present and must be ignored)

Links (both directions listed) with combined_score, threshold = 700.
"""
from __future__ import annotations

import os

import pandas as pd
import pyarrow.parquet as pq
import pytest

from string_network import load as L
from string_network.build import BuildOptions, build_network, write_output
from string_network.load import load_uniprot_idmapping

FIX = os.path.join(os.path.dirname(__file__), "fixtures")


def make_cfg():
    return {
        "taxon": 9606,
        "threshold": 700,
        "chunksize": 1000,
        "urls": {
            "links": os.path.join(FIX, "links.txt"),
            "uniprot_idmapping": os.path.join(FIX, "idmapping.dat"),
            "string_aliases": os.path.join(FIX, "aliases.txt"),
        },
        "paths": {"raw": FIX},
    }


def run(**over):
    cfg = make_cfg()
    opts = BuildOptions(
        threshold=over.pop("threshold", 700),
        threshold_when=over.pop("threshold_when", "after"),
        agg=over.pop("agg", "mean"),
        decimals=over.pop("decimals", 1),
        mapping_source="uniprot_idmapping",
        chunksize=1000,
    )
    df, stats = build_network(cfg, opts, verbose=False)
    return df, stats


def canon_pairs(df):
    return {
        tuple(sorted((a, b)))
        for a, b in zip(df["protein1"].tolist(), df["protein2"].tolist())
    }


def test_mapping_strips_isoform_and_filters_string_only():
    m = load_uniprot_idmapping(os.path.join(FIX, "idmapping.dat"))
    # Only STRING rows kept (the Ensembl 'QX' row is dropped).
    assert "QX" not in set(m["uniprot"])
    # Isoform suffix removed: P4-2 -> P4.
    assert "P4" in set(m["uniprot"])
    assert not m["uniprot"].str.contains(r"-\d+$").any()


def test_ab_and_ba_collapse_into_one_row():
    df, _ = run(agg="mean")
    pairs = list(zip(df["protein1"], df["protein2"]))
    # (P1,P2) must appear exactly once, and never as its reverse too.
    canon = [tuple(sorted(p)) for p in pairs]
    assert len(canon) == len(set(canon)), "duplicate undirected pairs"
    assert ("P1", "P2") in set(canon)


def test_self_loops_removed():
    df, stats = run(agg="mean")
    # ENSPA-ENSPA -> P1-P1 must be gone.
    assert not ((df["protein1"] == df["protein2"]).any())
    assert stats.counts["self_loops_removed"] == 1


def test_one_to_many_expansion_and_aggregation():
    # ENSPB -> {P2, P2B}: the A-B link must expand to (P1,P2) and (P1,P2B).
    df, _ = run(agg="mean", decimals=1)
    pairs = canon_pairs(df)
    assert ("P1", "P2") in pairs
    assert ("P1", "P2B") in pairs
    # (P1,P2) aggregates directed scores [900,900,755,755] -> mean 827.5
    score = df[(df.protein1.isin(["P1", "P2"])) & (df.protein2.isin(["P1", "P2"]))]
    val = float(score["combined_score"].iloc[0])
    assert abs(val - 827.5) < 1e-4


def test_aggregation_methods_differ():
    # (P1,P3) directed scores = [800,800,650,650].
    def pval(agg, when="after"):
        df, _ = run(agg=agg, decimals=None, threshold_when=when)
        row = df[((df.protein1 == "P1") & (df.protein2 == "P3")) |
                 ((df.protein1 == "P3") & (df.protein2 == "P1"))]
        return None if row.empty else float(row["combined_score"].iloc[0])

    assert abs(pval("mean") - 725.0) < 1e-4
    assert abs(pval("max") - 800.0) < 1e-4
    # min after aggregation = 650 -> below threshold -> pair dropped
    assert pval("min", when="after") is None


def test_threshold_applied_after_vs_before():
    # agg=min on (P1,P3) with scores [800,800,650,650]:
    #  - after:  min=650 -> dropped
    #  - before: 650 links filtered first -> [800,800] -> min=800 -> kept
    df_after, _ = run(agg="min", decimals=None, threshold_when="after")
    df_before, _ = run(agg="min", decimals=None, threshold_when="before")
    p_after = canon_pairs(df_after)
    p_before = canon_pairs(df_before)
    assert ("P1", "P3") not in p_after
    assert ("P1", "P3") in p_before


def test_output_dtypes_and_no_index(tmp_path):
    df, _ = run(agg="mean", decimals=1)
    assert str(df["protein1"].dtype) == "string"
    assert str(df["protein2"].dtype) == "string"
    assert str(df["combined_score"].dtype) == "float32"

    out = tmp_path / "out.parquet"
    write_output(df, str(out))
    schema = pq.ParquetFile(str(out)).schema_arrow
    names = [f.name for f in schema]
    assert not any(n.startswith("__index_level_") or n == "index" for n in names)


def test_sorted_descending_and_no_duplicates():
    df, _ = run(agg="mean", decimals=1)
    s = df["combined_score"].to_numpy()
    assert (s[:-1] >= s[1:]).all(), "not sorted descending"
    canon = [tuple(sorted(p)) for p in zip(df["protein1"], df["protein2"])]
    assert len(canon) == len(set(canon))


# ---------------------------------------------------------------------------
# score_mode = physical_else_scaled (the discovered reference rule)
# ---------------------------------------------------------------------------

def make_cfg_sc():
    return {
        "taxon": 9606,
        "threshold": 700,
        "chunksize": 1000,
        "urls": {
            "links": os.path.join(FIX, "links_sc.txt"),
            "links_physical": os.path.join(FIX, "links_phys_sc.txt"),
            "protein_info": os.path.join(FIX, "protein_info_sc.txt"),
            "uniprot_idmapping": os.path.join(FIX, "idmapping_sc.dat"),
            "string_aliases": os.path.join(FIX, "aliases_missing.txt"),
        },
        "paths": {"raw": FIX},
    }


def run_sc(**over):
    opts = BuildOptions(
        threshold=over.pop("threshold", 700),
        threshold_when="after",
        agg="max",
        decimals=over.pop("decimals", 1),
        mapping_source="uniprot_idmapping",
        score_mode=over.pop("score_mode", "physical_else_scaled"),
        scale_factor=over.pop("scale_factor", 0.9),
        chunksize=1000,
    )
    return build_network(make_cfg_sc(), opts, verbose=False)


def test_score_mode_physical_edge_uses_physical_score():
    # (Q1,Q2) is in the physical network -> uses physical score 985,
    # not the full score 990 nor 0.9*990.
    df, _ = run_sc()
    row = df[((df.protein1 == "Q1") & (df.protein2 == "Q2")) |
             ((df.protein1 == "Q2") & (df.protein2 == "Q1"))]
    assert abs(float(row["combined_score"].iloc[0]) - 985.0) < 1e-4


def test_score_mode_non_physical_edge_is_scaled_and_rounded():
    # (Q1,Q3) is NOT physical -> 0.9 * 901 = 810.9 (non-integer, 1 dp).
    df, _ = run_sc()
    row = df[((df.protein1 == "Q1") & (df.protein2 == "Q3")) |
             ((df.protein1 == "Q3") & (df.protein2 == "Q1"))]
    assert abs(float(row["combined_score"].iloc[0]) - 810.9) < 1e-4


def test_score_mode_threshold_applied_after_transform():
    # (Q2,Q3): full score 750 (>=700) but 0.9*750 = 675 (<700) -> dropped.
    df, _ = run_sc()
    pairs = {tuple(sorted(p)) for p in zip(df.protein1, df.protein2)}
    assert ("Q2", "Q3") not in pairs
    # raw mode would keep it (750 >= 700)
    df_raw, _ = run_sc(score_mode="raw")
    pairs_raw = {tuple(sorted(p)) for p in zip(df_raw.protein1, df_raw.protein2)}
    assert ("Q2", "Q3") in pairs_raw


def test_threshold_500_subset_equals_threshold_700():
    # Building at 500 and then filtering to >=700 must equal building at 700:
    # same pairs, same scores, same row count.
    df500, _ = run_sc(threshold=500)
    df700, _ = run_sc(threshold=700)
    sub = (df500[df500.combined_score >= 700]
           [["protein1", "protein2", "combined_score"]]
           .sort_values(["protein1", "protein2"]).reset_index(drop=True))
    exp = (df700[["protein1", "protein2", "combined_score"]]
           .sort_values(["protein1", "protein2"]).reset_index(drop=True))
    assert len(sub) == len(exp)
    assert sub.equals(exp)


def test_mapping_fallback_chain():
    # ENSP4 has no UniProt and no links; the fallback chain must still map it
    # to its preferred_name (G4), not drop it, when mapping_fallback is on.
    m = L.load_mapping_with_fallback(make_cfg_sc(), "uniprot_idmapping", 9606)
    d = dict(zip(m.ensp, m.uniprot))
    assert d["ENSP1"] == "Q1"          # UniProt accession wins
    assert d["ENSP4"] == "G4"          # falls back to preferred_name
