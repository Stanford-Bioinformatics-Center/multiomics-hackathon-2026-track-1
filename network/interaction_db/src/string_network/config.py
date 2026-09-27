"""Configuration loading and CLI-override helpers.

All defaults live in ``config.yaml``. Nothing else in the codebase hard-codes a
URL or a path; everything is read from the loaded config dict.
"""
from __future__ import annotations

import os
from typing import Any, Dict

import yaml

# Project root = two levels up from this file (src/string_network/config.py).
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def load_config(path: str = "config.yaml") -> Dict[str, Any]:
    """Load the YAML config and resolve all paths relative to the project root."""
    if not os.path.isabs(path):
        path = os.path.join(PROJECT_ROOT, path)
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    # Resolve paths relative to the project root so nothing is hard-coded.
    paths = cfg.setdefault("paths", {})
    for key, value in list(paths.items()):
        if value and not os.path.isabs(value):
            paths[key] = os.path.join(PROJECT_ROOT, value)

    cfg["_project_root"] = PROJECT_ROOT
    return cfg


def apply_overrides(cfg: Dict[str, Any], args: Any) -> Dict[str, Any]:
    """Override config values from parsed argparse args (only when not None)."""
    mapping = {
        "threshold": "threshold",
        "threshold_when": "threshold_when",
        "agg": "agg",
        "decimals": "decimals",
        "mapping_source": "mapping_source",
        "reviewed_only": "reviewed_only",
        "keep_list": "keep_list",
        "links_source": "links_source",
        "score_mode": "score_mode",
        "scale_factor": "scale_factor",
        "mapping_fallback": "mapping_fallback",
    }
    for arg_name, cfg_key in mapping.items():
        val = getattr(args, arg_name, None)
        if val is not None:
            cfg[cfg_key] = val

    # Comma-separated list override, e.g. --mapping-alias-sources UniProt_AC,Ensembl_UniProt
    alias = getattr(args, "mapping_alias_sources", None)
    if alias is not None:
        cfg["mapping_alias_sources"] = [s.strip() for s in alias.split(",") if s.strip()]
    return cfg


def resolve_output_path(cfg: Dict[str, Any], override: str = None) -> str:
    """Resolve the output parquet path.

    Precedence: an explicit *override* (e.g. --output) wins; otherwise the
    ``paths.output`` template is formatted with the (possibly overridden)
    threshold, so --threshold 500 writes ..._ge500.parquet while the default of
    700 keeps the historical ..._ge700.parquet path.
    """
    if override:
        return override
    template = cfg["paths"]["output"]
    return template.format(threshold=cfg.get("threshold", 700))
