"""ChEBI ontology + structure parsing.

chebi.obo carries names, synonyms, obsolete/replaced_by, is_a and relationships,
and property_value formula / charge / smiles, but NOT InChIKey. The 2025 ChEBI
download layout moved flat files to ``flat_files/``; InChI/InChIKey come from
``flat_files/structures.tsv.gz`` (columns: id, compound_id, status_id, molfile,
smiles, standard_inchi, standard_inchi_key, dimension, default_structure).
"""
from __future__ import annotations

import gzip
import re
from typing import Dict, List, Tuple

import pandas as pd

# The 2025 ChEBI obo encodes relationships with RO/chebi ontology ids rather than
# readable names. Map the ones we use to friendly names; readable names are kept
# too for backward compatibility with older obo files.
RO_MAP = {
    "RO:0018034": "is_conjugate_acid_of",
    "RO:0018033": "is_conjugate_base_of",
    "RO:0018036": "is_tautomer_of",
    "RO:0000087": "has_role",
}
REL_KEEP = {
    "is_conjugate_acid_of", "is_conjugate_base_of", "is_tautomer_of",
    "has_functional_parent", "has_major_microspecies_at_pH_7_3",
    "has_parent_hydride", "is_enantiomer_of",
}
# ChEBI's 2025 obo uses CURIEs like `chemrof:generalized_empirical_formula`,
# `chemrof:charge`, `chemrof:smiles_string` (older files used a URL path ending
# in /formula etc). Accept a ':' or '/' before the key.
PROP_RE = re.compile(r'property_value:\s+\S*?[:/](\w+)\s+"(.*?)"')
PROP_KEYS = {
    "formula": "formula",
    "generalized_empirical_formula": "formula",
    "empirical_formula": "formula",
    "charge": "charge",
    "smiles": "smiles",
    "smiles_string": "smiles",
}


def parse_obo(path: str):
    """Parse chebi.obo(.gz). Returns (names_df, relations_df, structures_df)."""
    op = gzip.open if path.endswith(".gz") else open
    names, relations, structures = [], [], []
    cur = None

    def flush(t):
        if not t or "id" not in t:
            return
        cid = t["id"]
        names.append({
            "chebi_id": cid,
            "name": t.get("name"),
            "synonyms": "|".join(t.get("synonyms", [])) or None,
            "is_obsolete": bool(t.get("is_obsolete")),
            "replaced_by": t.get("replaced_by"),
            "alt_ids": "|".join(t.get("alt_id", [])) or None,
        })
        structures.append({
            "chebi_id": cid,
            "formula": t.get("formula"),
            "charge": t.get("charge"),
            "smiles": t.get("smiles"),
        })
        for rel, obj in t.get("rels", []):
            relations.append({"subject": cid, "relation": rel, "object": obj})

    with op(path, "rt", encoding="utf-8", errors="replace") as fh:
        in_term = False
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                flush(cur)
                cur = {"synonyms": [], "alt_id": [], "rels": []}
                in_term = True
                continue
            if line.startswith("["):        # some other stanza
                flush(cur); cur = None; in_term = False
                continue
            if not in_term or cur is None or not line:
                continue
            if line.startswith("id: "):
                cur["id"] = line[4:].strip()
            elif line.startswith("name: "):
                cur["name"] = line[6:].strip()
            elif line.startswith("alt_id: "):
                cur["alt_id"].append(line[8:].strip())
            elif line.startswith("synonym: "):
                m = re.search(r'"(.*?)"', line)
                if m:
                    cur["synonyms"].append(m.group(1))
            elif line.startswith("is_obsolete: true"):
                cur["is_obsolete"] = True
            elif line.startswith("replaced_by: "):
                cur["replaced_by"] = line[13:].strip()
            elif line.startswith("consider: ") and not cur.get("replaced_by"):
                cur["replaced_by"] = line[10:].strip()
            elif line.startswith("is_a: "):
                obj = line[6:].split("!")[0].strip()
                cur["rels"].append(("is_a", obj))
            elif line.startswith("relationship: "):
                parts = line[len("relationship: "):].split("!")[0].split()
                if len(parts) >= 2:
                    rel = RO_MAP.get(parts[0], parts[0])
                    if rel in REL_KEEP:
                        cur["rels"].append((rel, parts[1]))
            elif line.startswith("property_value: "):
                m = PROP_RE.search(line)
                if m:
                    key = PROP_KEYS.get(m.group(1).lower())
                    if key and not cur.get(key):
                        cur[key] = m.group(2)
        flush(cur)

    names_df = pd.DataFrame(names)
    relations_df = pd.DataFrame(relations)
    structures_df = pd.DataFrame(structures)
    return names_df, relations_df, structures_df


def load_inchikeys(structures_tsv_gz: str) -> pd.DataFrame:
    """Parse flat_files/structures.tsv.gz -> chebi_id, inchikey, inchi.

    Keeps the default structure per compound. The molfile column contains
    multi-line quoted values, so the file is parsed in chunks with proper
    quoting and the molfile column is dropped immediately.
    """
    keep = ["compound_id", "standard_inchi", "standard_inchi_key", "default_structure"]
    out = []
    reader = pd.read_csv(
        structures_tsv_gz, sep="\t", dtype=str, quotechar='"',
        engine="c", chunksize=200000, on_bad_lines="skip",
    )
    for chunk in reader:
        chunk = chunk[[c for c in keep if c in chunk.columns]]
        chunk = chunk[chunk["default_structure"].astype(str).str.lower().isin(["true", "1", "t"])]
        out.append(chunk)
    df = pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=keep)
    df["chebi_id"] = "CHEBI:" + df["compound_id"].astype(str)
    df = df.rename(columns={"standard_inchi_key": "inchikey", "standard_inchi": "inchi"})
    df = df.dropna(subset=["inchikey"]).drop_duplicates("chebi_id")
    return df[["chebi_id", "inchikey", "inchi"]].reset_index(drop=True)


def alt_to_primary(names_df: pd.DataFrame) -> Dict[str, str]:
    """Map secondary (alt_id) and obsolete replaced_by ids to their primary id."""
    m: Dict[str, str] = {}
    for _, r in names_df.iterrows():
        if r.get("alt_ids"):
            for a in str(r["alt_ids"]).split("|"):
                if a:
                    m[a] = r["chebi_id"]
        if r.get("is_obsolete") and r.get("replaced_by"):
            m[r["chebi_id"]] = r["replaced_by"]
    return m


def run(cfg, args) -> int:
    print("[mnet.chebi] use the mnet-download step (download.run) which invokes "
          "these parsers.")
    return 0
