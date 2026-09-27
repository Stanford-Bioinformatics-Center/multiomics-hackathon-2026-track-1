"""The build pipeline: map -> dedup -> aggregate -> filter -> sort -> write.

Every filtering / merging step logs counts so nothing is dropped silently.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field, replace
from typing import Callable, Dict, List, Optional

import numpy as np
import pandas as pd

from . import load as L

AGG_FUNCS = {"mean", "max", "min", "median", "first"}


@dataclass
class BuildOptions:
    threshold: int = 700
    threshold_when: str = "after"          # "before" | "after"
    agg: str = "max"                        # mean | max | min | median | first
    decimals: Optional[int] = None
    mapping_source: str = "uniprot_idmapping"
    reviewed_only: bool = False
    keep_list: Optional[str] = None
    chunksize: int = 2_000_000
    links_source: str = "full"             # full | physical
    # score_mode:
    #   raw                  -> report the links_source combined score as-is
    #   physical_else_scaled -> physical score if the edge is in the physical
    #                           network, else scale_factor * full score
    score_mode: str = "raw"
    scale_factor: float = 0.9              # multiplier for non-physical edges
    mapping_fallback: bool = False         # UniProt -> preferred_name -> raw STRING id
    mapping_alias_sources: tuple = ("UniProt_AC", "Ensembl_UniProt")
    taxon: int = 9606

    @classmethod
    def from_config(cls, cfg: dict) -> "BuildOptions":
        return cls(
            threshold=int(cfg.get("threshold", 700)),
            threshold_when=str(cfg.get("threshold_when", "after")),
            agg=str(cfg.get("agg", "max")),
            decimals=cfg.get("decimals", None),
            mapping_source=str(cfg.get("mapping_source", "uniprot_idmapping")),
            reviewed_only=bool(cfg.get("reviewed_only", False)),
            keep_list=cfg.get("keep_list", None),
            chunksize=int(cfg.get("chunksize", 2_000_000)),
            links_source=str(cfg.get("links_source", "full")),
            score_mode=str(cfg.get("score_mode", "raw")),
            scale_factor=float(cfg.get("scale_factor", 0.9)),
            mapping_fallback=bool(cfg.get("mapping_fallback", False)),
            mapping_alias_sources=tuple(
                cfg.get("mapping_alias_sources", ("UniProt_AC", "Ensembl_UniProt"))
            ),
            taxon=int(cfg.get("taxon", 9606)),
        )


@dataclass
class BuildStats:
    """Counts logged at each step, returned alongside the result frame."""
    log: List[str] = field(default_factory=list)
    counts: Dict[str, int] = field(default_factory=dict)
    verbose: bool = True

    def add(self, key: str, value: int, msg: Optional[str] = None) -> None:
        self.counts[key] = value
        self.log.append(msg or f"{key}: {value:,}")
        if self.verbose:
            print(f"[build] {self.log[-1]}")


def _mapping_frame(cfg: dict, opts: BuildOptions, stats: BuildStats) -> pd.DataFrame:
    if opts.mapping_fallback:
        m = L.load_mapping_with_fallback(cfg, opts.mapping_source, opts.taxon,
                                         alias_sources=opts.mapping_alias_sources)
        stats.add("mapping_rows", len(m),
                  f"mapping rows ({opts.mapping_source} + preferred_name/raw fallback): {len(m):,}")
    else:
        m = L.load_mapping(cfg, opts.mapping_source)
        stats.add("mapping_rows", len(m), f"mapping rows ({opts.mapping_source}): {len(m):,}")
    stats.add("mapping_unique_ensp", m["ensp"].nunique())
    stats.add("mapping_unique_uniprot", m["uniprot"].nunique())

    # Explicit one-to-many diagnostics.
    ensp_fanout = m.groupby("ensp").size()
    uni_fanout = m.groupby("uniprot").size()
    stats.add("ensp_with_multi_uniprot", int((ensp_fanout > 1).sum()),
              f"ENSP mapping to >1 UniProt: {int((ensp_fanout > 1).sum()):,}")
    stats.add("uniprot_with_multi_ensp", int((uni_fanout > 1).sum()),
              f"UniProt mapped from >1 ENSP: {int((uni_fanout > 1).sum()):,}")
    return m


def links_url_key(links_source: str) -> str:
    """Map a links_source ('full'|'physical') to a config URL key."""
    return "links_physical" if links_source == "physical" else "links"


def links_path_for(cfg: dict, links_source: str = "full") -> str:
    """Resolve the path to the raw links file for the given source."""
    key = links_url_key(links_source)
    return os.path.join(cfg["paths"]["raw"], os.path.basename(cfg["urls"][key]))


def load_links_full(cfg: dict, chunksize: Optional[int] = None,
                    links_source: str = "full") -> pd.DataFrame:
    """Load the entire links file (no threshold). Useful for caching in tuning."""
    links_path = links_path_for(cfg, links_source)
    frames = list(L.iter_links(
        links_path, threshold=None,
        chunksize=chunksize or int(cfg.get("chunksize", 2_000_000)),
        apply_threshold=False,
    ))
    return pd.concat(frames, ignore_index=True)


def _load_links(cfg: dict, opts: BuildOptions, stats: BuildStats,
                links_cache: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    apply_before = opts.threshold_when == "before"
    if links_cache is not None:
        links = links_cache
        if apply_before and opts.threshold is not None:
            links = links[links["combined_score"] >= opts.threshold]
    else:
        links_path = links_path_for(cfg, opts.links_source)
        frames = list(L.iter_links(
            links_path,
            threshold=opts.threshold,
            chunksize=opts.chunksize,
            apply_threshold=apply_before,
        ))
        links = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
            columns=["protein1", "protein2", "combined_score"]
        )
    when = "before mapping (>= threshold)" if apply_before else "unfiltered"
    stats.add("links_loaded", len(links), f"directed links loaded ({when}): {len(links):,}")
    return links


def _map_ends(links: pd.DataFrame, m: pd.DataFrame, stats: BuildStats) -> pd.DataFrame:
    """Merge the mapping onto both ends; log link loss from unmapped ends."""
    m1 = m.rename(columns={"ensp": "protein1", "uniprot": "u1"})
    m2 = m.rename(columns={"ensp": "protein2", "uniprot": "u2"})

    merged = links.merge(m1, on="protein1", how="left")
    merged = merged.merge(m2, on="protein2", how="left")

    unmapped1 = merged["u1"].isna()
    unmapped2 = merged["u2"].isna()
    stats.add("links_unmapped_end1", int(unmapped1.sum()))
    stats.add("links_unmapped_end2", int(unmapped2.sum()))
    lost = (unmapped1 | unmapped2).sum()
    stats.add("links_lost_unmapped", int(lost),
              f"directed links lost (an end unmapped): {int(lost):,}")

    merged = merged.dropna(subset=["u1", "u2"])
    # One-to-many mapping naturally expands rows via the merge (Cartesian on ends).
    stats.add("links_after_mapping_expand", len(merged),
              f"directed links after mapping/expansion: {len(merged):,}")
    return merged[["u1", "u2", "combined_score"]]


def _dedup_aggregate(mapped: pd.DataFrame, opts: BuildOptions,
                     stats: BuildStats) -> pd.DataFrame:
    # Drop self-loops (same UniProt on both ends after mapping).
    self_mask = mapped["u1"] == mapped["u2"]
    stats.add("self_loops_removed", int(self_mask.sum()))
    mapped = mapped[~self_mask]

    # Canonical undirected key: sorted (a, b).
    u1 = mapped["u1"].to_numpy(dtype=object)
    u2 = mapped["u2"].to_numpy(dtype=object)
    a = np.where(u1 <= u2, u1, u2)
    b = np.where(u1 <= u2, u2, u1)
    canon = pd.DataFrame({
        "a": a,
        "b": b,
        "combined_score": mapped["combined_score"].to_numpy(),
    })
    stats.add("directed_rows_pre_agg", len(canon))

    agg = opts.agg
    if agg not in AGG_FUNCS:
        raise ValueError(f"Unknown agg: {agg!r}")

    grouped = canon.groupby(["a", "b"], sort=False)["combined_score"]
    if agg == "first":
        out = grouped.first()
    else:
        out = grouped.agg(agg)
    out = out.reset_index()
    stats.add("undirected_pairs", len(out),
              f"undirected pairs after {agg} aggregation: {len(out):,}")

    if opts.decimals is not None:
        out["combined_score"] = out["combined_score"].round(opts.decimals)
    else:
        # keep numeric; may be non-integer if agg == mean/median
        out["combined_score"] = out["combined_score"].astype("float64")

    return out.rename(columns={"a": "protein1", "b": "protein2"})


def _apply_optional_filters(df: pd.DataFrame, cfg: dict, opts: BuildOptions,
                            stats: BuildStats) -> pd.DataFrame:
    if opts.reviewed_only:
        reviewed = L.load_reviewed_set(cfg)
        keep = df["protein1"].isin(reviewed) & df["protein2"].isin(reviewed)
        stats.add("dropped_not_reviewed", int((~keep).sum()))
        df = df[keep]
        stats.add("pairs_after_reviewed_filter", len(df))
    if opts.keep_list:
        keep_set = L.load_keep_list(opts.keep_list)
        keep = df["protein1"].isin(keep_set) & df["protein2"].isin(keep_set)
        stats.add("dropped_not_in_keeplist", int((~keep).sum()))
        df = df[keep]
        stats.add("pairs_after_keeplist", len(df))
    return df


def _aggregate_physical(cfg: dict, opts: BuildOptions, m: pd.DataFrame,
                        stats: BuildStats,
                        physical_cache: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """Load + map + dedup the physical links network (max agg, no threshold).

    Returns canonical protein1/protein2/combined_score. Used by the score
    transform to decide which edges are 'physical'.
    """
    phys_opts = replace(opts, links_source="physical", threshold_when="after",
                        agg="max", decimals=None)
    links = _load_links(cfg, phys_opts, stats, links_cache=physical_cache)
    mapped = _map_ends(links, m, stats)
    pagg = _dedup_aggregate(mapped, phys_opts, stats)
    stats.add("physical_undirected_pairs", len(pagg))
    return pagg


def _apply_score_transform(cfg: dict, opts: BuildOptions, m: pd.DataFrame,
                           full_out: pd.DataFrame, stats: BuildStats,
                           physical_cache: Optional[pd.DataFrame] = None) -> pd.DataFrame:
    """score = physical score if edge in physical network, else scale_factor*full.

    ``full_out`` is the canonical full-network aggregation (no threshold yet).
    """
    phys = _aggregate_physical(cfg, opts, m, stats, physical_cache=physical_cache)
    merged = full_out.merge(
        phys.rename(columns={"combined_score": "phys_score"}),
        on=["protein1", "protein2"], how="outer",
    )
    full_s = merged["combined_score"].astype("float64")
    phys_s = merged["phys_score"].astype("float64")
    in_phys = phys_s.notna()

    score = np.where(in_phys, phys_s, opts.scale_factor * full_s)
    merged["combined_score"] = score
    decimals = 1 if opts.decimals is None else opts.decimals
    merged["combined_score"] = merged["combined_score"].round(decimals)

    stats.add("edges_from_physical", int(in_phys.sum()),
              f"edges taking the physical score: {int(in_phys.sum()):,}")
    stats.add("edges_scaled_0.9", int((~in_phys).sum()),
              f"edges scaled by {opts.scale_factor}: {int((~in_phys).sum()):,}")
    return merged[["protein1", "protein2", "combined_score"]]


def build_network(cfg: dict, opts: Optional[BuildOptions] = None,
                  links_cache: Optional[pd.DataFrame] = None,
                  mapping_cache: Optional[pd.DataFrame] = None,
                  physical_cache: Optional[pd.DataFrame] = None,
                  verbose: bool = True):
    """Run the full build. Returns (DataFrame, BuildStats).

    ``links_cache`` (full unfiltered links) and ``mapping_cache`` let callers
    such as the tuner avoid re-reading the large gzip files for every combo.
    """
    opts = opts or BuildOptions.from_config(cfg)
    stats = BuildStats(verbose=verbose)
    if verbose:
        print(f"[build] options: {opts}")

    if mapping_cache is not None:
        m = mapping_cache
        stats.add("mapping_rows", len(m))
    else:
        m = _mapping_frame(cfg, opts, stats)

    # When transforming scores we need the *untruncated* base scores, so the
    # threshold is only applied at the very end.
    base_opts = opts
    if opts.score_mode != "raw":
        base_opts = replace(opts, threshold_when="after")

    links = _load_links(cfg, base_opts, stats, links_cache=links_cache)
    mapped = _map_ends(links, m, stats)
    out = _dedup_aggregate(mapped, base_opts, stats)

    if opts.score_mode == "raw":
        pass
    elif opts.score_mode == "physical_else_scaled":
        out = _apply_score_transform(cfg, opts, m, out, stats,
                                     physical_cache=physical_cache)
    else:
        raise ValueError(f"Unknown score_mode: {opts.score_mode!r}")

    # Threshold is applied after aggregation and (if any) score transform.
    if opts.threshold_when == "after" or opts.score_mode != "raw":
        before = len(out)
        out = out[out["combined_score"] >= opts.threshold]
        stats.add("dropped_below_threshold_after_agg", before - len(out))
        stats.add("pairs_after_threshold", len(out))

    out = _apply_optional_filters(out, cfg, opts, stats)

    # Sort by score descending and finalise dtypes.
    out = out.sort_values("combined_score", ascending=False, kind="mergesort")
    out = out.reset_index(drop=True)
    out["protein1"] = out["protein1"].astype("string")
    out["protein2"] = out["protein2"].astype("string")
    out["combined_score"] = out["combined_score"].astype("float32")

    nodes = pd.unique(pd.concat([out["protein1"], out["protein2"]], ignore_index=True))
    stats.add("final_edges", len(out))
    stats.add("final_nodes", len(nodes))
    if len(out) and verbose:
        s = out["combined_score"]
        print(f"[build] score min/median/max: {s.min():.1f} / "
              f"{s.median():.1f} / {s.max():.1f}")
    return out, stats


def write_output(df: pd.DataFrame, path: str) -> None:
    """Write parquet with no pandas index via pyarrow."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    df.to_parquet(path, engine="pyarrow", index=False)
    print(f"[build] wrote {path} ({len(df):,} rows)")
