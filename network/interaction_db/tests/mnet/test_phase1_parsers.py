"""Hermetic unit tests for Phase 1 pure parsing logic (no network, no big files)."""
from __future__ import annotations

from mnet import chebi as chebi_mod
from mnet import rhea as rhea_mod
from mnet.download import _norm_name


def test_classify_compound_type():
    f = rhea_mod.classify_compound_type
    assert f("C10H16N5O13P3") == "small molecule"   # ATP
    assert f(None) == "other"
    assert f("") == "other"
    # R-group placeholder -> generic
    assert f("C5H9NO2R") == "generic"
    # polymer repeat notation -> polymer
    assert f("(C6H10O5)n") == "polymer"


def test_chebi_prop_regex_new_and_old_formats():
    # 2025 CURIE style
    m = chebi_mod.PROP_RE.search(
        'property_value: chemrof:generalized_empirical_formula "C6H12O6" xsd:string')
    assert m and chebi_mod.PROP_KEYS[m.group(1).lower()] == "formula"
    assert m.group(2) == "C6H12O6"
    m = chebi_mod.PROP_RE.search('property_value: chemrof:charge "-4" xsd:integer')
    assert m and chebi_mod.PROP_KEYS[m.group(1).lower()] == "charge" and m.group(2) == "-4"
    # legacy URL style still handled
    m = chebi_mod.PROP_RE.search(
        'property_value: http://purl.obolibrary.org/obo/chebi/formula "H2O" xsd:string')
    assert m and chebi_mod.PROP_KEYS[m.group(1).lower()] == "formula"


def test_norm_name_strips_charge_and_zwitterion():
    assert _norm_name("ATP(4-)") == "atp"
    assert _norm_name("S-adenosyl-L-methionine zwitterion") == "s-adenosyl-l-methionine"
    assert _norm_name(float("nan")) == ""


def test_short_chebi_and_side():
    assert rhea_mod._short_chebi("http://purl.obolibrary.org/obo/CHEBI_16459") == "CHEBI:16459"
    assert rhea_mod._side("http://rdf.rhea-db.org/10000_R") == "R"
    assert rhea_mod._master_id("http://rdf.rhea-db.org/10000") == 10000
