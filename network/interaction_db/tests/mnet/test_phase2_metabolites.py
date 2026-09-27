"""Phase 2 tests: charge normalization, combined rows, obsolete, manual mappings."""
from __future__ import annotations

import os

import pandas as pd
import pytest

from mnet import metabolites as M
from mnet.metabolites import Ctx
from string_network.config import load_config

ROOT = os.path.join(os.path.dirname(__file__), "..", "..")
CONFIG = os.path.join(ROOT, "config.yaml")
INTERIM = os.path.join(ROOT, "data", "interim")


def _bare_ctx(**attrs) -> Ctx:
    c = object.__new__(Ctx)
    c.ph7 = {}
    c.adj = {}
    c.all_rhea_chebi = set()
    c.human_rhea_chebi = set()
    c.human_rxn_count = {}
    c.obsolete = {}
    c.name_index = {}
    c.__dict__.update(attrs)
    return c


# ---------------- hermetic ----------------

def test_normalize_ph7_primary():
    c = _bare_ctx(ph7={"CHEBI:15366": "CHEBI:30089"})
    assert c.normalize("CHEBI:15366") == ("CHEBI:30089", "ph7_3_mapping")


def test_normalize_conjugate_fallback_prefers_rhea():
    c = _bare_ctx(adj={"CHEBI:1": {"CHEBI:2", "CHEBI:3"}}, all_rhea_chebi={"CHEBI:3"})
    cid, method = c.normalize("CHEBI:1")
    assert cid == "CHEBI:3" and method == "conjugate_tautomer"


def test_normalize_none_when_no_mapping():
    c = _bare_ctx()
    assert c.normalize("CHEBI:999") == ("CHEBI:999", "none")


def test_deobsolete():
    c = _bare_ctx(obsolete={"CHEBI:old": "CHEBI:new"})
    assert c.deobsolete("CHEBI:old") == ("CHEBI:new", True)
    assert c.deobsolete("CHEBI:keep") == ("CHEBI:keep", False)


def test_refmet_node_fallback():
    assert M._refmet_node("RM0001", "x") == "REFMET:RM0001"
    assert M._refmet_node(float("nan"), "Leucine/Isoleucine") == "MEAS:Leucine_Isoleucine"
    assert M._refmet_node(None, "a b") == "MEAS:a_b"


def test_load_manual_reads_filled_skips_blank(tmp_path):
    p = tmp_path / "manual.csv"
    p.write_text("metabolite,refmet_id,suggested_chebi,note\n"
                 "Foo,RM1,CHEBI:123,ok\n"
                 "Bar,RM2,,todo\n")
    m = M._load_manual(str(p))
    assert m == {"RM1": "CHEBI:123"}


# ---------------- data-dependent (skip if interim caches absent) ----------------

_need = ["chebi_names", "chebi_structures", "chebi_relations",
         "rhea_participants", "rhea_enzymes"]
_have_interim = all(os.path.exists(os.path.join(INTERIM, f"{n}.parquet")) for n in _need)
skip_no_data = pytest.mark.skipif(not _have_interim, reason="interim caches not built")


@pytest.fixture(scope="module")
def real_ctx():
    return Ctx(load_config(CONFIG))


@skip_no_data
def test_charge_normalization_real(real_ctx):
    # expected ChEBI ids verified against Rhea's chebi_pH7_3 mapping
    assert real_ctx.normalize("CHEBI:15366")[0] == "CHEBI:30089"   # acetic acid -> acetate
    assert real_ctx.normalize("CHEBI:422")[0] == "CHEBI:16651"     # L-lactic -> (S)-lactate
    assert real_ctx.normalize("CHEBI:30915")[0] == "CHEBI:16810"   # 2-oxoglutaric -> 2-oxoglutarate(2-)


@skip_no_data
def test_combined_row_two_chebi():
    o = pd.read_parquet(os.path.join(INTERIM, "metabolites_nonlipid.parquet"))
    r = o[o.metabolite == "Leucine/Isoleucine"]
    assert len(r) == 1
    assert len(str(r.iloc[0]["chebi_rhea"]).split(";")) == 2


@skip_no_data
def test_manual_mapping_applied():
    # write a manual override for a known unmapped metabolite, run once, assert
    # applied, then restore all outputs from backups (no second pipeline run).
    manual = os.path.join(ROOT, "curation", "manual_metabolite_mappings.csv")
    outputs = [
        os.path.join(INTERIM, "metabolites_nonlipid.parquet"),
        os.path.join(ROOT, "reports", "mnet_phase2_metabolites.md"),
        os.path.join(ROOT, "reports", "mnet_phase2_mapping_review.csv"),
    ]
    backups = {manual: open(manual, "rb").read()}
    for p in outputs:
        backups[p] = open(p, "rb").read() if os.path.exists(p) else None
    try:
        o = pd.read_parquet(outputs[0])
        unm = o[(o.mapping_method == "unmapped") & o.refmet_id.notna()]
        if unm.empty:
            pytest.skip("no unmapped row with refmet to test")
        target = unm.iloc[0]
        with open(manual, "w") as fh:
            fh.write("metabolite,refmet_id,suggested_chebi,note\n")
            fh.write(f"{target['metabolite']},{target['refmet_id']},CHEBI:16810,test\n")
        M.run(load_config(CONFIG), None)
        o2 = pd.read_parquet(outputs[0])
        row = o2[o2.refmet_id == target["refmet_id"]].iloc[0]
        assert row["mapping_method"] == "manual"
        assert "CHEBI:16810" in str(row["chebi_rhea"])
    finally:
        for p, data in backups.items():
            if data is not None:
                with open(p, "wb") as fh:
                    fh.write(data)
