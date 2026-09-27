"""Load the existing STRING protein-protein network as ppi edges.

Will read the STRING parquet at ``string_network.config.resolve_output_path(cfg)``
(unless ``mnet.ppi_network`` overrides it) and emit ppi edges unchanged: same
pairs, same scores.

Phase 0: scaffold only.
"""
from __future__ import annotations

import os


def ppi_path(cfg) -> str:
    """Resolve the STRING PPI parquet used by the metabolite network.

    Uses mnet.ppi_network if set, else the string_network output template
    formatted with mnet.ppi_threshold (default: the string_network threshold).
    """
    mcfg = cfg.get("mnet", {}) or {}
    p = mcfg.get("ppi_network")
    if p:
        return p if os.path.isabs(p) else os.path.join(cfg["_project_root"], p)
    thr = mcfg.get("ppi_threshold", cfg.get("threshold", 700))
    return cfg["paths"]["output"].format(threshold=thr)


def load_ppi(cfg):
    import pandas as pd
    return pd.read_parquet(ppi_path(cfg))


def ppi_nodes(cfg) -> set:
    df = load_ppi(cfg)
    return set(df["protein1"]) | set(df["protein2"])


def run(cfg, args) -> int:
    print("[mnet.network] Phase 0 scaffold - not implemented yet.")
    return 0
