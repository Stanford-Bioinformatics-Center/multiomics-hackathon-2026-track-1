"""Phase 3 lipid tests (small, per SPEED MODE)."""
from __future__ import annotations

import os

import pandas as pd
import pytest

from mnet import lipids as LP

INTERIM = os.path.join(os.path.dirname(__file__), "..", "..", "data", "interim")


def test_parse_and_chain_helpers():
    cls, level, tc, tdb, tox, chains = LP._parse_lipid("PC 34:1")
    assert cls == "PC" and level == "species" and tc == 34 and tdb == 1
    # carbon-from-name + chain class
    assert LP._carbon_from_name("Octadecenoic acid") == 18
    assert LP._carbon_from_name("FA 16:0") == 16
    assert LP._chain_class(4) == "short"
    assert LP._chain_class(10) == "medium"
    assert LP._chain_class(18) == "long"
    assert LP._chain_class(26) == "very-long"


def test_canon_matches_swisslipids_and_pygoslin_forms():
    # SwissLipids "PC(34:1)" and pygoslin "PC 34:1" must canonicalise equal
    assert LP._canon("PC(34:1)") == LP._canon("PC 34:1") == "PC34:1"


@pytest.mark.skipif(
    not os.path.exists(os.path.join(INTERIM, "edges_lipid.parquet")),
    reason="phase 3 outputs not built")
def test_lipid_is_a_edges():
    e = pd.read_parquet(os.path.join(INTERIM, "edges_lipid.parquet"))
    assert (e["edge_type"] == "lipid_is_a").all()
    assert e["node2"].str.startswith("LIPIDCLASS:").all()
    assert (e["node1"] != e["node2"]).all()          # no self-loops
