"""argparse entry point: ``python -m mnet <command>``.

Phase 0 scaffold: every subcommand is wired to its module's ``run`` function,
which currently reports that the step is not implemented yet. No data is
processed. Configuration is read from the shared ``config.yaml`` via
``string_network.config.load_config``; the mnet-specific settings live under the
top-level ``mnet:`` key.
"""
from __future__ import annotations

import argparse
import sys
from typing import Any, Dict

from string_network.config import load_config

from . import (assemble, chebi, download, lipids, metabolites, network,
               proteins, ptm, qc, rhea)

# subcommand -> module (each module exposes run(cfg, args))
COMMANDS = {
    "download": download,
    "metabolites": metabolites,
    "lipids": lipids,
    "rhea": rhea,
    "proteins": proteins,
    "network": network,
    "assemble": assemble,
    "qc": qc,
    "ptm": ptm,
}


def mnet_cfg(cfg: Dict[str, Any]) -> Dict[str, Any]:
    """Return the ``mnet:`` sub-config (empty dict if absent)."""
    return cfg.get("mnet", {}) or {}


def build_parser() -> argparse.ArgumentParser:
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--config", default="config.yaml")

    ap = argparse.ArgumentParser(prog="mnet", parents=[common],
                                 description="Metabolite extension of the STRING network.")
    sub = ap.add_subparsers(dest="command", required=True)
    for name in COMMANDS:
        sub.add_parser(name, parents=[common], help=f"{name} step")
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_config(args.config)
    module = COMMANDS[args.command]
    return int(module.run(cfg, args) or 0)


if __name__ == "__main__":
    sys.exit(main())
