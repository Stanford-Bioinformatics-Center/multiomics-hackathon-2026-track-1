"""Phase 6 PTM tests (max 3, per SPEED MODE)."""
from __future__ import annotations

import os

import pandas as pd
import pytest

from mnet import ptm

OUT = os.path.join(os.path.dirname(__file__), "..", "..", "data", "output")


def test_uniprot_feature_parse(tmp_path):
    p = tmp_path / "u.tsv"
    seq = "MKKSPTYAAANQTGGG"          # S4 T6 Y7 ...; N12 sequon N-Q-T
    p.write_text(
        "Entry\tGene Names (primary)\tSequence\tModified residue\tGlycosylation\n"
        "TESTP\tG1\t" + seq + "\t"
        'MOD_RES 4; /note="Phosphoserine; by CDK1 and GSK3B"; /evidence="ECO:0000269|PubMed:1"; '
        'MOD_RES 7; /note="Phosphotyrosine"; /evidence="ECO:0000250|UniProtKB:X"\t'
        'CARBOHYD 12; /note="N-linked (GlcNAc...) asparagine"; /evidence="ECO:0000269|PubMed:2"\n')
    up = ptm._parse_uniprot(str(p))
    rec = up["TESTP"]
    assert rec["gene"] == "G1"
    ph = {(pos, res): (evi, kin) for pos, res, note, evi, kin in rec["phospho"]}
    assert ph[(4, "S")] == ("experimental", ["CDK1", "GSK3B"])
    assert ph[(7, "Y")][0] == "predicted/similarity"     # ECO:0000250
    assert rec["glyco"][0][2] == "N-linked"


def test_motrpac_parse_and_isoform_mapping():
    # single / multi / isoform feature_id parsing
    assert ptm.MOTRPAC_RE.match("P46020_S758s")
    m = ptm.MOTRPAC_RE.match("A0FGR8-2_S665sS676sT677t")
    assert m and m.group(1) == "A0FGR8-2"
    assert ptm.SITE_RE.findall(m.group(2)) == [("S", "665"), ("S", "676"), ("T", "677")]
    assert ptm.MOTRPAC_RE.match("NOTAFEATURE") is None
    # isoform window -> canonical position (unique 15-mer, centre residue at index 7)
    window = "ABCDEFGHIJKLMNO"
    canon = "XXABCDEFGHIJKLMNOZZ"
    pos, status = ptm._map_isoform(window, canon)
    assert status == "isoform_mapped" and canon[pos - 1] == "H"


@pytest.mark.skipif(not os.path.exists(os.path.join(OUT, "proteins_ptm.csv")),
                    reason="phase 6 outputs not built")
def test_per_protein_aggregation_zeros():
    df = pd.read_csv(os.path.join(OUT, "proteins_ptm.csv"))
    # counts are integers with 0 (not blank) for proteins with no sites
    for c in ["n_phosphosites", "n_glycosites", "n_substrate_sites", "n_motrpac_sites"]:
        assert df[c].notna().all()
        assert (df[c] >= 0).all()
    nosite = df[df.n_phosphosites == 0]
    assert len(nosite) > 0 and (nosite["n_phosphosites"] == 0).all()
