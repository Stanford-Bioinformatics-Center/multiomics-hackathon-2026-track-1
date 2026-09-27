"""Phase 4 Rhea-layer tests (data-dependent; skip if interim not built)."""
from __future__ import annotations

import os

import pandas as pd
import pytest

INTERIM = os.path.join(os.path.dirname(__file__), "..", "..", "data", "interim")
_have = os.path.exists(os.path.join(INTERIM, "edges_rhea.parquet"))
skip = pytest.mark.skipif(not _have, reason="phase 4 outputs not built")

GLUCOSE = "CHEBI:4167"           # D-glucopyranose
G6P = "CHEBI:58225"              # D-glucose 6-phosphate(2-)
ATP = "CHEBI:30616"


@skip
def test_no_metabolite_metabolite_edges():
    e = pd.read_parquet(os.path.join(INTERIM, "edges_rhea.parquet"))
    # every Rhea-layer edge is protein -> metabolite/class
    assert (e["node1_type"] == "protein").all()
    assert e["node2_type"].isin(["metabolite", "lipid_class"]).all()


@skip
def test_hexokinase_atp_dropped_but_kept_for_synthesis():
    e = pd.read_parquet(os.path.join(INTERIM, "edges_rhea.parquet"))
    atp = e[e["node2"] == ATP]
    # kinase reaction RHEA:17825 (glucose + ATP -> G6P + ADP): ATP is a cofactor -> dropped
    assert not atp["evidence"].str.contains("RHEA:17825").any()
    # but ATP is still kept for reactions where it is a real substrate (synthesis/ligases)
    assert len(atp) > 0


@skip
def test_crosswalk_lands_on_existing_string_nodes():
    cross = pd.read_parquet(os.path.join(INTERIM, "protein_crosswalk.parquet"))
    from mnet import network as net
    from string_network.config import load_config
    string_nodes = net.ppi_nodes(load_config(
        os.path.join(os.path.dirname(__file__), "..", "..", "config.yaml")))
    # accessions resolved via ENSP / gene must land on a real STRING node
    via = cross[cross["method"].isin(["via_ensp", "via_gene"]) & cross["in_string"]]
    assert len(via) >= 1
    assert via["node_id"].isin(string_nodes).all()
    # and the vast majority of enzymes are on existing STRING nodes
    assert cross["in_string"].mean() > 0.9
