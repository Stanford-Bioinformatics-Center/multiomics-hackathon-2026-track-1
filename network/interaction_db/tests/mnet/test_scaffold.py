"""Phase 0 scaffold tests for the mnet package.

These must not require the heavy mnet dependencies (rdflib, pygoslin, ...) or any
downloaded data - they only check that the package and CLI are wired up and that
the shared config exposes the mnet settings.
"""
from __future__ import annotations

import os

import pytest

import mnet
from mnet import cli
from string_network.config import load_config

CONFIG = os.path.join(os.path.dirname(__file__), "..", "..", "config.yaml")

EXPECTED_COMMANDS = {
    "download", "metabolites", "lipids", "rhea", "proteins",
    "network", "assemble", "qc",
}


def test_package_imports():
    assert mnet.__version__


def test_cli_parser_has_all_commands():
    assert EXPECTED_COMMANDS.issubset(set(cli.COMMANDS))
    parser = cli.build_parser()
    assert parser is not None


# Only the commands whose run() is still a harmless print-only stub. The data
# steps (download, metabolites, lipids, rhea, assemble, qc) are implemented and
# do real work / network I/O, so they are not exercised here.
STUB_COMMANDS = {"proteins", "network"}


def test_stub_commands_run(capsys):
    cfg = load_config(CONFIG)

    class _Args:
        config = CONFIG

    for name in STUB_COMMANDS:
        rc = cli.COMMANDS[name].run(cfg, _Args())
        assert rc == 0, f"{name}.run should return 0"


def test_mnet_config_section_present():
    cfg = load_config(CONFIG)
    m = cfg.get("mnet")
    assert m is not None, "config.yaml must have an mnet: section"
    for key in ("metabolites", "interim", "output_dir", "human_taxon",
                "include_unmeasured", "keep_currency_in_catalysis", "lipid_is_a"):
        assert key in m, f"mnet.{key} missing"


def test_new_url_keys_present():
    cfg = load_config(CONFIG)
    urls = cfg["urls"]
    for key in ("rhea_rhea2uniprot", "rhea_directions", "rhea_rdf",
                "chebi_obo", "swisslipids_lipids"):
        assert key in urls, f"urls.{key} missing"
    # existing keys must remain
    for key in ("links", "links_physical", "uniprot_idmapping", "string_aliases"):
        assert key in urls, f"existing urls.{key} must be unchanged"
