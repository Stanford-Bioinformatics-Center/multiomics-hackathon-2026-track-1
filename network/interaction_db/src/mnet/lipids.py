"""Lipids: parse, map to classes, build lipid -> class (lipid_is_a) edges.

Reuses the Phase-2 metabolites.Ctx for the heavy ChEBI / Rhea lookups. pygoslin
parses the complex-lipid and acylcarnitine shorthand; common-name fatty acids are
handled by carbon count -> chain-length ChEBI class. SwissLipids is matched by a
fast string-canonicalised join on the Abbreviation column (not pygoslin on 779k
rows).
"""
from __future__ import annotations

import os
import re
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from . import metabolites as MB

LIPID_SUPERCLASSES = {"Glycerolipids", "Glycerophospholipids", "Sphingolipids"}

# Rhea generic ChEBI classes actually present in human Rhea (anion / pH 7.3 forms).
# Each FA chain class = its specific chain-length anion (if in human Rhea) PLUS the
# general fatty-acid anion (shared across all FA classes).
FA_ANION_GENERAL = "CHEBI:28868"        # 'fatty acid anion'
FA_CHAIN_CHEBI = {
    "short": ["CHEBI:28868"],                          # no short-chain anion in human Rhea
    "medium": ["CHEBI:59558", "CHEBI:28868"],          # medium-chain fatty acid anion
    "long": ["CHEBI:57560", "CHEBI:28868"],            # long-chain fatty acid anion
    "very-long": ["CHEBI:58950", "CHEBI:28868"],       # very long-chain fatty acid anion
}
CAR_CHEBI = ("CHEBI:75659", "O-acyl-L-carnitine")

# amino acids etc. that are NOT fatty acids and must not enter the lipid scope
NON_FATTY_ACID = {"aminobutyric acid"}

# carbon-count prefixes for common fatty-acid names
C_PREFIX = {
    "propan": 3, "propen": 3, "butan": 4, "buten": 4, "butyr": 4, "pentan": 5,
    "penten": 5, "hexan": 6, "hexen": 6, "heptan": 7, "hepten": 7, "octan": 8,
    "octen": 8, "nonan": 9, "nonen": 9, "decan": 10, "decen": 10, "dodec": 12,
    "tridec": 13, "tetradec": 14, "pentadec": 15, "hexadec": 16, "heptadec": 17,
    "octadec": 18, "nonadec": 19, "eicos": 20, "icos": 20, "heneicos": 21,
    "docos": 22, "tricos": 23, "tetracos": 24, "pentacos": 25, "hexacos": 26,
}
GENERAL_STOP_NAMES = {
    "lipid", "fatty acid", "fatty acid anion", "carboxylic acid",
    "glycerophospholipid", "glycerolipid", "phospholipid", "sphingolipid",
    "glyceride", "fatty acyl-coa", "membrane lipid", "polar lipid",
}

# Substrings identifying a class's GENERIC ChEBI (R-group forms) among human-Rhea
# participants - deliberately require "acyl"/generic wording so specific acyl
# species (e.g. 1-hexadecanoyl-...) are NOT swept in. Capped per class below.
CLASS_NAME_HINTS = {
    "PC": ["diacyl-sn-glycero-3-phosphocholine"],
    "PE": ["diacyl-sn-glycero-3-phosphoethanolamine"],
    "PS": ["diacyl-sn-glycero-3-phospho-l-serine"],
    "PI": ["diacyl-sn-glycero-3-phospho-1d-myo-inositol",
           "diacyl-sn-glycero-3-phospho-(1d-myo-inositol)"],
    "PG": ["diacyl-sn-glycero-3-phospho-(1'-sn-glycerol)"],
    "PA": ["diacyl-sn-glycero-3-phosphate", "1,2-diacyl-sn-glycero-3-phosphate"],
    "TG": ["triacyl-sn-glycerol", "tri(acyl)-sn-glycerol"],
    "DG": ["diacyl-sn-glycerol", "1,2-diacyl-sn-glycerol"],
    "MG": ["monoacyl-sn-glycerol", "monoacylglycerol"],
    "SM": ["acylsphingosine-1-phosphocholine", "sphingomyelin"],
    "Cer": ["n-acylsphingosine", "n-acylsphinganine", "n-(acyl)sphing"],
    "LPC": ["1-acyl-sn-glycero-3-phosphocholine", "2-acyl-sn-glycero-3-phosphocholine",
            "monoacyl-sn-glycero-3-phosphocholine", "1-o-acyl-sn-glycero-3-phosphocholine"],
    "LPE": ["1-acyl-sn-glycero-3-phosphoethanolamine",
            "2-acyl-sn-glycero-3-phosphoethanolamine",
            "monoacyl-sn-glycero-3-phosphoethanolamine"],
    "SPB": ["sphinganine", "sphingosine", "sphing-4-enine"],
    "SPBP": ["sphinganine 1-phosphate", "sphingosine 1-phosphate",
             "sphing-4-enine 1-phosphate"],
}
# ether (plasmanyl/plasmenyl = alkyl/alkenyl) generics for the *-ether classes
ETHER_HINTS = {
    "PC-ether": ["alkyl-sn-glycero-3-phosphocholine", "alkenyl-sn-glycero-3-phosphocholine"],
    "PE-ether": ["alkyl-sn-glycero-3-phosphoethanolamine", "alkenyl-sn-glycero-3-phosphoethanolamine"],
    "LPC-ether": ["alkyl-sn-glycero-3-phosphocholine", "alkenyl-sn-glycero-3-phosphocholine"],
    "LPE-ether": ["alkyl-sn-glycero-3-phosphoethanolamine", "alkenyl-sn-glycero-3-phosphoethanolamine"],
}
MAX_CLASS_CHEBI = 12


def _canon(s: str) -> str:
    """Canonicalise a lipid shorthand for matching (drop spaces/parens)."""
    if not isinstance(s, str):
        return ""
    return s.replace("(", "").replace(")", "").replace(" ", "").upper()


def _chain_class(c: int) -> str:
    if c < 6:
        return "short"
    if c <= 12:
        return "medium"
    if c <= 22:
        return "long"
    return "very-long"


def _carbon_from_name(name: str) -> Optional[int]:
    n = name.lower()
    m = re.search(r"\bfa\s+(\d+):", n) or re.search(r"\bcar\s+(\d+):", n)
    if m:
        return int(m.group(1))
    best = None
    for pref, c in C_PREFIX.items():
        if pref in n:
            best = c if best is None else max(best, c)
    return best


def _parse_lipid(name: str):
    """pygoslin parse -> (lipid_class, level, total_c, total_db, total_ox, chains)."""
    from pygoslin.parser.Parser import LipidParser
    from pygoslin.domain.LipidLevel import LipidLevel
    global _PARSER
    try:
        _PARSER
    except NameError:
        _PARSER = LipidParser()
    la = _PARSER.parse(name)
    lip = la.lipid
    cls = la.get_lipid_string(LipidLevel.CLASS)
    info = lip.info
    level = info.level.name.lower() if info and info.level else "undefined"
    def _fa_db(fa):
        v = getattr(fa, "double_bonds", None)
        try:
            v = v() if callable(v) else v
            return int(v)
        except Exception:
            return "?"
    chains = []
    try:
        for fa in lip.fa_list:
            chains.append(f"{fa.num_carbon}:{_fa_db(fa)}")
    except Exception:
        pass
    def _asint(v):
        try:
            v = v() if callable(v) else v
            return int(v)
        except Exception:
            return None
    tc = _asint(getattr(info, "num_carbon", None)) if info else None
    tdb = _asint(getattr(info, "double_bonds", None)) if info else None
    if tdb is None and info is not None:
        tdb = _asint(getattr(info, "get_double_bonds", None))
    tox = _asint(getattr(info, "num_oxygens", None)) if info else None
    return cls, level, tc, tdb, tox, chains


def _build_swisslipids_index(sl: pd.DataFrame) -> Dict[str, Dict]:
    """canon(Abbreviation) -> {slm, chebi, sl_class, level} (prefer rows with CHEBI)."""
    idx: Dict[str, Dict] = {}
    ab = sl["Abbreviation*"].fillna("")
    for slm, a, chebi, lvl, cls in zip(sl["Lipid ID"], ab, sl["CHEBI"],
                                       sl["Level"], sl["Lipid class*"]):
        key = _canon(a)
        if not key:
            continue
        rec = {"slm": slm, "chebi": chebi if isinstance(chebi, str) and chebi else None,
               "level": lvl, "sl_class": cls}
        if key not in idx or (rec["chebi"] and not idx[key]["chebi"]):
            idx[key] = rec
    return idx


def _isa_maps(cfg) -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]]]:
    rel = MB._read(cfg, "chebi_relations")
    isa = rel[rel["relation"] == "is_a"]
    up: Dict[str, Set[str]] = {}
    down: Dict[str, Set[str]] = {}
    for s, o in zip(isa["subject"], isa["object"]):
        up.setdefault(s, set()).add(o)
        down.setdefault(o, set()).add(s)
    return up, down


def run(cfg, args) -> int:
    ctx = MB.Ctx(cfg)
    interim = MB._interim(cfg)
    curation_dir = os.path.join(cfg["_project_root"], "curation")

    # participant compound_type map
    parts = MB._read(cfg, "rhea_participants")
    ctype = (parts.sort_values("compound_type")
             .drop_duplicates("chebi_id").set_index("chebi_id")["compound_type"].to_dict())

    csv = pd.read_csv((cfg.get("mnet", {}) or {}).get("metabolites"), dtype=str) \
        if os.path.isabs((cfg.get("mnet", {}) or {}).get("metabolites", "")) \
        else pd.read_csv(os.path.join(cfg["_project_root"],
                         (cfg.get("mnet", {}) or {}).get("metabolites")), dtype=str)
    phase2 = MB._read(cfg, "metabolites_nonlipid")
    p2_by_ref = phase2.set_index("refmet_id", drop=False)

    sl = MB._read(cfg, "swisslipids")
    sl_idx = _build_swisslipids_index(sl)

    # ---- scope ----
    lip_rows = csv[csv["super_class"].isin(LIPID_SUPERCLASSES)].copy()
    lip_rows["kind"] = "complex_lipid"
    # ALL Fatty Acyls (incl. specific FAs like Margaric/Pentadecylic) so each gets
    # an is_a edge to its chain-length class, plus any CAR rows.
    fa_car = csv[(csv["super_class"] == "Fatty Acyls") |
                 csv["metabolite"].str.startswith("CAR ", na=False)].copy()
    # drop non-fatty-acids that slipped in (e.g. Aminobutyric acid) - kept in the
    # manual mappings file for the reviewer to decide.
    fa_car = fa_car[~fa_car["metabolite"].str.strip().str.lower().isin(NON_FATTY_ACID)]
    fa_car["kind"] = "fa_car"
    scope = pd.concat([lip_rows, fa_car], ignore_index=True).drop_duplicates("metabolite")

    parse_fail = []
    out_rows = []
    for _, r in scope.iterrows():
        name = str(r["metabolite"]).strip()
        refmet = r.get("refmet_id")
        node_fallback = MB._refmet_node(refmet, name)
        alt = None
        lipid_class = level = None
        tc = tdb = tox = None
        chains = []
        exact_chebi = None
        mapping_method = "class_only"
        confidence = "medium"

        # either/or ether (e.g. "PE P-40:6 or PE O-40:7")
        cand_names = [name]
        if " or " in name:
            cand_names = [x.strip() for x in name.split(" or ")]
            alt = ";".join(cand_names)

        parsed = None
        for cn in cand_names:
            try:
                parsed = _parse_lipid(cn)
                break
            except Exception:
                continue
        if parsed:
            lipid_class, level, tc, tdb, tox, chains = parsed
            if len(cand_names) > 1:
                lipid_class = f"{lipid_class}-ether"
        else:
            parse_fail.append(name)

        # reuse Phase-2 normalized ChEBI (complex lipids with a ChEBI, and specific
        # fatty acids) as the exact node. class_level fatty acids keep a REFMET node
        # (their class ChEBI is only an attribute; they get a class is_a edge).
        if refmet in set(p2_by_ref.index):
            p2 = p2_by_ref.loc[refmet]
            if isinstance(p2, pd.DataFrame):
                p2 = p2.iloc[0]
            cr = p2.get("chebi_rhea")
            if pd.notna(cr) and str(cr) and ";" not in str(cr) and \
               p2.get("mapping_method") != "class_level":
                exact_chebi = str(cr)
                mapping_method = "phase2_chebi"

        # SwissLipids match (species / molecular levels)
        slm = None
        if parsed:
            from pygoslin.domain.LipidLevel import LipidLevel
            for lvl in (LipidLevel.SPECIES, LipidLevel.MOLECULAR_SPECIES):
                try:
                    key = _canon(_PARSER.parse(cand_names[0]).get_lipid_string(lvl))
                except Exception:
                    continue
                if key in sl_idx:
                    rec = sl_idx[key]
                    slm = rec["slm"]
                    if rec["chebi"]:
                        cid = rec["chebi"] if rec["chebi"].startswith("CHEBI:") else f"CHEBI:{rec['chebi']}"
                        if cid in ctx.human_rhea_chebi and exact_chebi is None:
                            exact_chebi = cid
                            mapping_method = "swisslipids_chebi"
                    break

        # ---- fatty acids / acylcarnitines ----
        is_car = name.upper().startswith("CAR ")
        is_fa = (r["kind"] == "fa_car") and not is_car
        class_node = None
        class_chebis: List[str] = []
        if is_car:
            c = _carbon_from_name(name)
            lipid_class = "acylcarnitine"
            class_node = "LIPIDCLASS:acylcarnitine"
            # specific O-palmitoyl-L-carnitine for CAR 16:0
            if c == 16 and ":0" in name and "CHEBI:17490" in ctx.human_rhea_chebi:
                exact_chebi = "CHEBI:17490"; mapping_method = "specific_chebi"
        elif is_fa:
            c = _carbon_from_name(name)
            chain = _chain_class(c) if c else "long"
            lipid_class = f"FA-{chain}-chain"
            class_node = f"LIPIDCLASS:{lipid_class}"
        else:
            # complex lipid -> class node from pygoslin class
            if lipid_class:
                class_node = f"LIPIDCLASS:{lipid_class}"
        # class_chebi_ids for every class are filled from the class table below

        # node id
        if exact_chebi:
            node_id = f"CHEBI:{exact_chebi.split(':')[-1]}"
            confidence = "high"
        else:
            node_id = node_fallback

        out_rows.append({
            "node_id": node_id, "metabolite": name, "refmet_id": refmet,
            "lipid_class": lipid_class, "lipid_level": level,
            "chains": ";".join(chains) if chains else None,
            "total_c": tc, "total_db": tdb, "total_oxygens": tox,
            "swisslipids_id": slm, "exact_chebi": exact_chebi,
            "class_node_id": class_node,
            "class_chebi_ids": ";".join(class_chebis) if class_chebis else None,
            "alt_forms": alt, "mapping_method": mapping_method, "confidence": confidence,
        })

    out = pd.DataFrame(out_rows)

    # distinct human reactions per chebi (for accurate class reaction counts)
    hp = parts[parts["master_rhea_id"].isin(
        set(MB._read(cfg, "rhea_enzymes")
            .assign(h=lambda d: d["uniprot"].isin(ctx._human_set()))
            .query("h")["MASTER_ID"].dropna().astype(int)))]
    chebi2rx = hp.groupby("chebi_id")["master_rhea_id"].apply(set).to_dict()

    def ndistinct(ids):
        s = set()
        for c in ids:
            s |= chebi2rx.get(c, set())
        return len(s)

    # ---- class -> ChEBI table (all class types) ----
    class_chebi_rows, claim = _class_to_chebi(cfg, out, sl, ctx, ctype, curation_dir)
    cls_map = {c: ids for c, ids in class_chebi_rows}
    # fill complex-lipid class_chebi_ids from the table
    def fill_cls(row):
        existing = row["class_chebi_ids"]
        if isinstance(existing, str) and existing:
            return existing
        node = row["class_node_id"]
        if isinstance(node, str) and node.startswith("LIPIDCLASS:"):
            cls = node.split(":", 1)[1]
            return ";".join(cls_map.get(cls, [])) or None
        return None
    out["class_chebi_ids"] = out.apply(fill_cls, axis=1)

    # ---- edges: lipid -> class ----
    edges = out[out["class_node_id"].notna()][["node_id", "class_node_id"]].copy()
    edges = edges[edges["node_id"] != edges["class_node_id"]].drop_duplicates()
    edges.columns = ["node1", "node2"]
    edges["node1_type"] = "metabolite"
    edges["node2_type"] = "lipid_class"
    edges["edge_type"] = "lipid_is_a"

    # ---- lipid_classes table (distinct human reactions) ----
    lc_rows = []
    for node, sub in out[out.class_node_id.notna()].groupby("class_node_id"):
        cls = node.split(":", 1)[1]
        ids = cls_map.get(cls, [])
        lc_rows.append({
            "class_node_id": node, "class_chebi_ids": ";".join(ids) or None,
            "n_lipids": sub["node_id"].nunique(),
            "n_human_rhea_reactions": ndistinct(ids),
            "n_shared_generics": sum(1 for c in ids if claim.get(c, 0) > 1),
        })
    lipid_classes = pd.DataFrame(lc_rows)

    # ---- write ----
    out.to_parquet(os.path.join(interim, "lipids.parquet"), index=False)
    lipid_classes.to_parquet(os.path.join(interim, "lipid_classes.parquet"), index=False)
    edges.to_parquet(os.path.join(interim, "edges_lipid.parquet"), index=False)
    _report(cfg, out, lipid_classes, edges, parse_fail, class_chebi_rows)
    _prefix_report(cfg, lipid_classes, claim, ctx)
    print(f"[mnet.lipids] lipids={len(out):,} classes={len(lipid_classes):,} "
          f"edges={len(edges):,} parse_fail={len(parse_fail)}")
    return 0


def _prefix_report(cfg, lipid_classes, claim, ctx):
    """Phase 4 pre-fix report: FA before/after, ether mapping, shared generics."""
    L = ["# mnet Phase 4 pre-fixes", ""]
    L.append("## 1. Non-fatty-acids removed from lipid scope")
    L.append("- Removed: **Aminobutyric acid** (amino acid, not a fatty acid); kept in "
             "curation/manual_metabolite_mappings.csv. The 4 hydroxy fatty acids "
             "(Hydroxyoctanoic/decanoic/dodecanoic/tetradecanoic) are genuine FAs and kept.")
    L.append("")
    L.append("## 2. FA chain-length classes: before -> after (distinct human reactions)")
    L.append("")
    before = {"LIPIDCLASS:FA-short-chain": 60, "LIPIDCLASS:FA-medium-chain": 60,
              "LIPIDCLASS:FA-long-chain": 8, "LIPIDCLASS:FA-very-long-chain": 60}
    L.append("| class | chebi_ids | before | after |")
    L.append("|---|---|---|---|")
    for node in ["LIPIDCLASS:FA-short-chain", "LIPIDCLASS:FA-medium-chain",
                 "LIPIDCLASS:FA-long-chain", "LIPIDCLASS:FA-very-long-chain"]:
        r = lipid_classes[lipid_classes.class_node_id == node]
        if len(r):
            rr = r.iloc[0]
            L.append(f"| {node} | {rr['class_chebi_ids']} | {before.get(node,'-')} | "
                     f"{rr['n_human_rhea_reactions']} |")
    L.append("- specific chain-length anions now used (medium CHEBI:59558, long "
             "CHEBI:57560, very-long CHEBI:58950) PLUS the shared general fatty-acid "
             "anion CHEBI:28868; short-chain has no specific anion in human Rhea.")
    L.append("")
    L.append("## 3. Ether classes")
    L.append("")
    for node in ["LIPIDCLASS:PC-ether", "LIPIDCLASS:PE-ether",
                 "LIPIDCLASS:LPC-ether", "LIPIDCLASS:LPE-ether"]:
        r = lipid_classes[lipid_classes.class_node_id == node]
        if len(r):
            rr = r.iloc[0]
            L.append(f"- {node}: {rr['n_human_rhea_reactions']} human reactions, "
                     f"chebis={rr['class_chebi_ids']}")
    L.append("(alkyl/alkenyl generics used where present; see how_found/confidence "
             "in curation/lipid_class_to_chebi.csv)")
    L.append("")
    L.append("## 4. Shared generic ChEBI (claimed by >1 class node)")
    L.append("")
    shared = {c: n for c, n in claim.items() if n > 1}
    L.append(f"- shared generic ChEBIs: **{len(shared)}**")
    for c, n in sorted(shared.items(), key=lambda x: -x[1]):
        L.append(f"  - {c} ({ctx.chebi_name.get(c,'?')}) claimed by {n} class nodes")
    total_shared_edges = int(lipid_classes["n_shared_generics"].sum())
    L.append(f"- total (class, shared-generic) claims across classes: "
             f"**{total_shared_edges}** (each generates Rhea catalysis edges to EACH "
             f"claiming class node in Phase 4, marked shared_generic=True)")
    L.append("")
    with open(os.path.join(cfg["paths"]["reports"], "mnet_phase4_prefixes.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def _class_to_chebi(cfg, out, sl, ctx, ctype, curation_dir):
    """Build/reuse class -> generic ChEBI table.

    Returns (list[(class, [chebi,...])], claim_counter). Handles complex classes
    (generic name hints; ether classes prefer alkyl/alkenyl generics, else diacyl
    fallback), FA chain-length classes and acylcarnitine. Marks shared_generic
    (a ChEBI claimed by >1 class node).
    """
    from collections import Counter
    path = os.path.join(curation_dir, "lipid_class_to_chebi.csv")

    def _claims(result):
        claim = Counter()
        for _cls, ids in result.items():
            for c in ids:
                claim[c] += 1
        return claim

    if os.path.exists(path) and os.path.getsize(path) > 60:  # user-edited version wins
        df = pd.read_csv(path, dtype=str)
        res = {cls: [c for c in sub["chebi_id"].dropna().tolist() if c]
               for cls, sub in df.groupby("class")}
        return list(res.items()), _claims(res)

    stop = set()
    for nm in GENERAL_STOP_NAMES:
        stop |= ctx.name_index.get(nm, set())

    sl_cls = sl[sl["Level"] == "Class"].copy()
    sl_cls["ab"] = sl_cls["Abbreviation*"].fillna("").str.replace(
        r"\(.*", "", regex=True).str.strip().str.upper()
    sl_class_chebi: Dict[str, Set[str]] = {}
    for ab, chebi in zip(sl_cls["ab"], sl_cls["CHEBI"]):
        if isinstance(chebi, str) and chebi:
            cid = chebi if chebi.startswith("CHEBI:") else f"CHEBI:{chebi}"
            sl_class_chebi.setdefault(ab, set()).add(cid)

    gen_names = {c: ctx.chebi_name.get(c, "").lower()
                 for c in ctx.human_rhea_chebi
                 if ctype.get(c) in ("generic", "small molecule") and ctx.chebi_name.get(c)}

    def _hint_hits(hints):
        hit = [c for c, nm in gen_names.items() if any(h in nm for h in hints)]
        return sorted(hit, key=lambda c: ctx.human_rxn_count.get(c, 0), reverse=True)[:MAX_CLASS_CHEBI]

    all_classes = sorted({n.split(":", 1)[1] for n in
                          out["class_node_id"].dropna().unique()})

    result: Dict[str, List[str]] = {}
    meta: Dict[str, Dict[str, Tuple[str, str]]] = {}   # cls -> {chebi: (how_found, confidence)}
    for cls in all_classes:
        found: Dict[str, Tuple[str, str]] = {}
        if cls.startswith("FA-") and cls.endswith("-chain"):
            chain = cls[len("FA-"):-len("-chain")]
            for c in FA_CHAIN_CHEBI.get(chain, [FA_ANION_GENERAL]):
                if c in ctx.human_rhea_chebi:
                    how = "fa_general" if c == FA_ANION_GENERAL else "fa_chain_specific"
                    found[c] = (how, "medium")
            if not found and FA_ANION_GENERAL in ctx.human_rhea_chebi:
                found[FA_ANION_GENERAL] = ("fa_general", "low")
        elif cls == "acylcarnitine":
            if CAR_CHEBI[0] in ctx.human_rhea_chebi:
                found[CAR_CHEBI[0]] = ("carnitine_generic", "medium")
        else:
            base = cls.replace("-ether", "")
            # ether classes: prefer alkyl/alkenyl generics
            if cls.endswith("-ether") and cls in ETHER_HINTS:
                for c in _hint_hits(ETHER_HINTS[cls]):
                    found[c] = ("ether_generic", "medium")
            if not found:
                howf = "diacyl_fallback" if cls.endswith("-ether") else "name_generic"
                conf = "low" if cls.endswith("-ether") else "high"
                for c in _hint_hits(CLASS_NAME_HINTS.get(base, [])):
                    found[c] = (howf, conf)
            for c in sl_class_chebi.get(base.upper(), set()):
                if c in ctx.human_rhea_chebi:
                    found.setdefault(c, ("swisslipids_class", "medium"))
            keep = {c for c in found if c not in stop} or set(found)
            found = {c: found[c] for c in keep}
        result[cls] = sorted(found)
        meta[cls] = found

    claim = _claims(result)
    rows = []
    for cls in all_classes:
        for c in result[cls]:
            how, conf = meta[cls][c]
            rows.append({"class": cls, "chebi_id": c,
                         "chebi_name": ctx.chebi_name.get(c, "?"),
                         "n_human_reactions": ctx.human_rxn_count.get(c, 0),
                         "how_found": how, "confidence": conf,
                         "shared_generic": claim[c] > 1})
    pd.DataFrame(rows, columns=["class", "chebi_id", "chebi_name", "n_human_reactions",
                                "how_found", "confidence", "shared_generic"]
                 ).to_csv(path, index=False)
    return list(result.items()), claim


def _report(cfg, out, lipid_classes, edges, parse_fail, class_chebi_rows):
    L = ["# mnet Phase 3 - lipids", ""]
    n = len(out)
    L.append(f"- rows in scope: {n} | parse success: {n-len(parse_fail)} "
             f"({100*(n-len(parse_fail))/n:.0f}%) | parse fail: {len(parse_fail)}")
    L.append(f"- SwissLipids hits: {out['swisslipids_id'].notna().sum()} | "
             f"exact_chebi assigned: {out['exact_chebi'].notna().sum()}")
    L.append(f"- lipid_is_a edges: {len(edges)} | class nodes: {len(lipid_classes)}")
    L.append("")
    L.append("## Counts by level")
    L.append(f"{out['lipid_level'].value_counts(dropna=False).to_dict()}")
    L.append("")
    L.append("## Counts by class (top 20)")
    L.append("")
    vc = out["lipid_class"].value_counts().head(20)
    L.append("| class | n |")
    L.append("|---|---|")
    for k, v in vc.items():
        L.append(f"| {k} | {v} |")
    L.append("")
    L.append("## Class -> ChEBI (with human Rhea reactions)")
    L.append("")
    L.append("| class_node | n_chebi | n_human_rhea_reactions |")
    L.append("|---|---|---|")
    for _, r in lipid_classes.sort_values("n_human_rhea_reactions", ascending=False).iterrows():
        nc = len(str(r["class_chebi_ids"]).split(";")) if r["class_chebi_ids"] else 0
        L.append(f"| {r['class_node_id']} | {nc} | {r['n_human_rhea_reactions']} |")
    L.append("")
    zero = lipid_classes[lipid_classes["n_human_rhea_reactions"] == 0]["class_node_id"].tolist()
    L.append(f"## Classes with 0 human Rhea reactions ({len(zero)})")
    L.append(f"{sorted(zero)}")
    L.append("")
    L.append(f"## Parse failures ({len(parse_fail)})")
    L.append(f"{sorted(parse_fail)[:40]}")
    L.append("")
    with open(os.path.join(cfg["paths"]["reports"], "mnet_phase3_lipids.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
