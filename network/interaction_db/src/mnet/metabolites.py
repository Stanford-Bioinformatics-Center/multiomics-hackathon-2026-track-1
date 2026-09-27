"""Map measured (non-lipid) metabolites to Rhea-ready ChEBI ids.

See Phase 2 brief. Produces data/interim/metabolites_nonlipid.parquet, a review
CSV and reports/mnet_phase2_metabolites.md. Complex lipids that already carry a
ChEBI id are normalised too (is_complex_lipid=True) for Phase 3 to reuse.
"""
from __future__ import annotations

import os
import re
from collections import deque
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from . import network as net

LIPID_SUPERCLASSES = {"Glycerolipids", "Glycerophospholipids", "Sphingolipids"}
FA_NOPOS_RE = re.compile(
    r"(?:^|[^0-9ZE])(?:mono|di|tri|tetra|penta|hexa)?en?oic acid$", re.I)
LOCANT_RE = re.compile(r"\d+[ZE]")


# --------------------------------------------------------------------------- IO
def _interim(cfg) -> str:
    p = (cfg.get("mnet", {}) or {}).get("interim", "data/interim")
    return p if os.path.isabs(p) else os.path.join(cfg["_project_root"], p)


def _read(cfg, name) -> pd.DataFrame:
    return pd.read_parquet(os.path.join(_interim(cfg), f"{name}.parquet"))


class Ctx:
    """Precomputed lookups shared across rows."""

    def __init__(self, cfg):
        self.cfg = cfg
        names = _read(cfg, "chebi_names")
        struct = _read(cfg, "chebi_structures")
        rel = _read(cfg, "chebi_relations")
        parts = _read(cfg, "rhea_participants")
        enz = _read(cfg, "rhea_enzymes")

        # obsolete -> replaced_by
        obs = names[(names["is_obsolete"]) & names["replaced_by"].notna()]
        self.obsolete = dict(zip(obs["chebi_id"], obs["replaced_by"]))

        # name / synonym -> set(chebi)
        self.name_index: Dict[str, Set[str]] = {}
        for _, r in names.iterrows():
            if r["name"]:
                self.name_index.setdefault(str(r["name"]).lower(), set()).add(r["chebi_id"])
            if r.get("synonyms"):
                for syn in str(r["synonyms"]).split("|"):
                    if syn:
                        self.name_index.setdefault(syn.lower(), set()).add(r["chebi_id"])
        self.chebi_name = dict(zip(names["chebi_id"], names["name"]))

        # inchikey lookups
        st = struct.dropna(subset=["inchikey"])
        self.ik_full: Dict[str, str] = dict(zip(st["inchikey"], st["chebi_id"]))
        self.ik_block1: Dict[str, Set[str]] = {}
        for ik, cid in zip(st["inchikey"], st["chebi_id"]):
            self.ik_block1.setdefault(ik[:14], set()).add(cid)

        # pH 7.3 mapping (numeric -> numeric)
        ph = pd.read_csv(os.path.join(cfg["paths"]["raw"], "chebi_pH7_3_mapping.tsv"),
                         sep="\t", dtype=str)
        ph.columns = [c.strip().upper() for c in ph.columns]
        self.ph7 = {f"CHEBI:{a}": f"CHEBI:{b}"
                    for a, b in zip(ph["CHEBI"], ph["CHEBI_PH7_3"])}

        # conjugate / tautomer adjacency (undirected)
        self.adj: Dict[str, Set[str]] = {}
        rr = rel[rel["relation"].isin(
            ["is_conjugate_acid_of", "is_conjugate_base_of", "is_tautomer_of"])]
        for s, o in zip(rr["subject"], rr["object"]):
            self.adj.setdefault(s, set()).add(o)
            self.adj.setdefault(o, set()).add(s)

        # is_a children (for class-level candidate isomers)
        isa = rel[rel["relation"] == "is_a"]
        self.children: Dict[str, List[str]] = {}
        for s, o in zip(isa["subject"], isa["object"]):
            self.children.setdefault(o, []).append(s)

        # human reactions + participant chebi counts
        human = self._human_set()
        human_reacs = set(enz[enz["uniprot"].isin(human)]["MASTER_ID"].dropna().astype(int))
        self.all_rhea_chebi: Set[str] = set(parts["chebi_id"])
        hp = parts[parts["master_rhea_id"].isin(human_reacs)]
        self.human_rhea_chebi: Set[str] = set(hp["chebi_id"])
        self.human_rxn_count: Dict[str, int] = (
            hp.groupby("chebi_id")["master_rhea_id"].nunique().to_dict())

        # currency chebi ids
        cur = pd.read_csv(os.path.join(cfg["_project_root"], "curation",
                                       "currency_metabolites.csv"), dtype=str)
        self.currency: Set[str] = set(cur["chebi_id"].dropna()) - {""}

    def _human_set(self) -> Set[str]:
        raw = self.cfg["paths"]["raw"]
        ac = set(pd.read_csv(os.path.join(raw, "HUMAN_9606_idmapping.dat.gz"),
                             sep="\t", header=None, usecols=[0], dtype=str)[0])
        with open(os.path.join(raw, f"reviewed_{self.cfg['taxon']}.list")) as fh:
            ac |= {ln.strip() for ln in fh if ln.strip()}
        return ac

    # -- helpers --
    def deobsolete(self, cid: str) -> Tuple[str, bool]:
        if cid in self.obsolete:
            return self.obsolete[cid], True
        return cid, False

    def normalize(self, cid: str) -> Tuple[str, str]:
        """Return (normalized_chebi, method). pH7.3 primary, conj/tautomer fallback."""
        if cid in self.ph7:
            return self.ph7[cid], "ph7_3_mapping"
        # BFS up to 3 hops; prefer a node that appears in Rhea participants
        seen = {cid}
        q = deque([(cid, 0)])
        cands = []
        while q:
            node, d = q.popleft()
            if d >= 3:
                continue
            for nb in self.adj.get(node, ()):
                if nb not in seen:
                    seen.add(nb)
                    cands.append(nb)
                    q.append((nb, d + 1))
        rhea_cands = [c for c in cands if c in self.all_rhea_chebi]
        if rhea_cands:
            best = max(rhea_cands, key=lambda c: (c in self.human_rhea_chebi,
                                                  self.human_rxn_count.get(c, 0)))
            return best, "conjugate_tautomer"
        return cid, "none"

    def rank(self, cid: str) -> Tuple[int, int]:
        return (int(cid in self.human_rhea_chebi), self.human_rxn_count.get(cid, 0))

    def best_rhea_form(self, cid: str) -> str:
        """Upgrade a (possibly general) ChEBI to the human-Rhea form actually used.

        Considers the id itself, its pH 7.3 form, its is_a children and their
        pH 7.3 forms (e.g. general 'leucine' -> L-leucine -> L-leucine zwitterion),
        and conjugate/tautomer neighbours. Prefers a form in human Rhea.
        """
        cands = {cid}
        if cid in self.ph7:
            cands.add(self.ph7[cid])
        for ch in self.children.get(cid, []):
            cands.add(ch)
            if ch in self.ph7:
                cands.add(self.ph7[ch])
        cands |= self.adj.get(cid, set())
        in_rhea = [c for c in cands if c in self.human_rhea_chebi] or \
                  [c for c in cands if c in self.all_rhea_chebi]
        if in_rhea:
            return max(in_rhea, key=self.rank)
        return self.ph7.get(cid, cid)

    def best_rhea_forms(self, cid: str, max_forms: int = 4):
        """All human-Rhea forms of *cid* (pH7.3 + is_a children + conjugates).

        Returns the list of human-Rhea ChEBI ids; for an ambiguous parent this may
        be several (e.g. phosphoglycerate -> 2- and 3-; 2-hydroxyglutarate -> L/D).
        Returns [] if none, or if the parent is too ambiguous (> max_forms) so the
        caller can fall back (used to avoid exploding e.g. a fatty-acid class into
        all positional isomers).
        """
        cands = {cid}
        if cid in self.ph7:
            cands.add(self.ph7[cid])
        for ch in self.children.get(cid, []):
            cands.add(ch)
            if ch in self.ph7:
                cands.add(self.ph7[ch])
        cands |= self.adj.get(cid, set())
        forms = sorted(c for c in cands if c in self.human_rhea_chebi)
        if not forms or len(forms) > max_forms:
            return []
        return forms


def _pick(ctx: Ctx, norm_ids: List[str], log: List[str], label: str) -> str:
    """Pick the best normalized candidate; log ties."""
    if len(norm_ids) == 1:
        return norm_ids[0]
    ranked = sorted(set(norm_ids), key=ctx.rank, reverse=True)
    top = ctx.rank(ranked[0])
    tied = [c for c in ranked if ctx.rank(c) == top]
    if len(tied) > 1:
        log.append(f"[tie] {label}: {tied} all rank {top}; chose {sorted(tied)[0]}")
        return sorted(tied)[0]
    return ranked[0]


def _resolve_no_chebi(ctx: Ctx, row) -> Tuple[Optional[str], str, str, Optional[str]]:
    """Return (chebi, mapping_method, confidence, candidate_isomers)."""
    ik = row.get("inchi_key")
    name = str(row.get("metabolite") or "").strip()
    # a) full inchikey
    if isinstance(ik, str) and ik in ctx.ik_full:
        return ctx.ik_full[ik], "inchikey", "high", None
    # b) inchikey block1
    if isinstance(ik, str) and len(ik) >= 14:
        blk = ctx.ik_block1.get(ik[:14])
        if blk:
            if len(blk) == 1:
                return next(iter(blk)), "inchikey_block1", "medium", None
            pref = [c for c in blk if c in ctx.all_rhea_chebi]
            chosen = max(pref or blk, key=ctx.rank)
            return chosen, "inchikey_block1", "medium", ";".join(sorted(blk))
    # d) fatty acid without double-bond position -> class level
    nm = name.lower()
    if FA_NOPOS_RE.search(nm) and not LOCANT_RE.search(nm):
        hit = ctx.name_index.get(nm)
        if hit:
            cid = max(hit, key=ctx.rank)
            kids = ctx.children.get(cid, [])
            return cid, "class_level", "low", ";".join(sorted(kids)[:20]) or None
    # c) exact name / synonym
    hit = ctx.name_index.get(nm)
    if hit:
        cid = max(hit, key=ctx.rank)
        return cid, "name", "medium", None
    return None, "unmapped", "none", None


def _refmet_node(refmet, name: str) -> str:
    """REFMET:<id> when a refmet id exists, else a stable name-based MEAS id."""
    if refmet is not None and pd.notna(refmet) and str(refmet).strip().lower() not in ("", "nan"):
        return f"REFMET:{refmet}"
    slug = re.sub(r"[^A-Za-z0-9]+", "_", str(name)).strip("_")
    return f"MEAS:{slug}"


def _load_manual(path: str) -> Dict[str, str]:
    """refmet_id -> suggested_chebi for filled-in rows."""
    if not os.path.exists(path):
        return {}
    df = pd.read_csv(path, dtype=str)
    if "suggested_chebi" not in df.columns or "refmet_id" not in df.columns:
        return {}
    df = df[df["suggested_chebi"].notna() & (df["suggested_chebi"].str.strip() != "")]
    return dict(zip(df["refmet_id"], df["suggested_chebi"].str.strip()))


def run(cfg, args) -> int:
    ctx = Ctx(cfg)
    src = (cfg.get("mnet", {}) or {}).get("metabolites", "data/input/metabolite_chebi_ids.csv")
    src = src if os.path.isabs(src) else os.path.join(cfg["_project_root"], src)
    df = pd.read_csv(src, dtype=str)

    curation_dir = os.path.join(cfg["_project_root"], "curation")
    manual = _load_manual(os.path.join(curation_dir, "manual_metabolite_mappings.csv"))

    log: List[str] = []
    obsolete_replacements = 0

    # scope: non-lipid OR complex-lipid-with-chebi
    in_scope = (~df["super_class"].isin(LIPID_SUPERCLASSES)) | \
               (df["super_class"].isin(LIPID_SUPERCLASSES) & df["chebi_id"].notna())
    work = df[in_scope].copy()

    out_rows = []
    for _, row in work.iterrows():
        refmet = row.get("refmet_id")
        name = str(row.get("metabolite") or "").strip()
        is_lip = row["super_class"] in LIPID_SUPERCLASSES
        norm_method = "none"
        candidate_isomers = None

        # ---- manual override wins ----
        if refmet in manual:
            base_ids = [c.strip() for c in manual[refmet].split(";") if c.strip()]
            norm_ids, methods = [], []
            for b in base_ids:
                n, m = ctx.normalize(b)
                norm_ids.append(n); methods.append(m)
            norm_method = methods[0] if methods else "none"
            mapping_method, confidence = "manual", "medium"
            chosen = norm_ids
        else:
            # ---- combined two-metabolite row (e.g. Leucine/Isoleucine) ----
            combined = None
            if "/" in name and not is_lip:
                parts = [p.strip() for p in name.split("/")]
                resolved = []
                for p in parts:
                    hit = ctx.name_index.get(p.lower())
                    if hit:
                        resolved.append(max(hit, key=ctx.rank))
                if len(resolved) == len(parts) and len(parts) >= 2:
                    combined = resolved
            if combined is not None:
                norm_ids = []
                for b in combined:
                    n = ctx.best_rhea_form(b)     # upgrade to the human-Rhea form
                    norm_ids.append(n)
                norm_method = "ph7_3_mapping"
                mapping_method, confidence = "exact", "high"
                chosen = norm_ids
            else:
                # ---- candidates from chebi_id + chebi_all ----
                cands = []
                for col in ("chebi_id", "chebi_all"):
                    v = row.get(col)
                    if isinstance(v, str) and v.strip():
                        cands += [c.strip() for c in v.split(";") if c.strip()]
                cands = [c for c in dict.fromkeys(cands) if c.startswith("CHEBI:")]
                # de-obsolete
                deob = []
                for c in cands:
                    nc, changed = ctx.deobsolete(c)
                    if changed:
                        obsolete_replacements += 1
                        log.append(f"[obsolete] {c} -> {nc} ({name})")
                    deob.append(nc)
                cands = list(dict.fromkeys(deob))

                if cands:
                    norm_ids, methods, changed_any = [], [], False
                    for c in cands:
                        n, m = ctx.normalize(c)
                        norm_ids.append(n); methods.append(m)
                        if n != c:
                            changed_any = True
                    norm_method = ("ph7_3_mapping" if "ph7_3_mapping" in methods else
                                   ("conjugate_tautomer" if "conjugate_tautomer" in methods else "none"))
                    best = _pick(ctx, norm_ids, log, name)
                    chosen = [best]
                    mapping_method = "charge_normalized" if best != cands[0] and changed_any else "exact"
                    confidence = "high"
                else:
                    cid, mapping_method, confidence, candidate_isomers = _resolve_no_chebi(ctx, row)
                    if cid is None:
                        chosen = []
                    else:
                        n, norm_method = ctx.normalize(cid)
                        chosen = [n]

        # upgrade any measured id not used in human Rhea to the form(s) actually
        # used (pH7.3 / conjugate / is_a children). Ambiguous parents (e.g.
        # phosphoglycerate -> 2-/3-, 2-hydroxyglutarate -> L/D) attach all forms.
        # class_level fatty acids expand to too many isomers -> best_rhea_forms
        # returns [] and they keep their class-level id (they get a class is_a edge).
        upgraded, changed = [], False
        for c in chosen:
            if c in ctx.human_rhea_chebi:
                upgraded.append(c)
            else:
                forms = ctx.best_rhea_forms(c)
                if forms:
                    upgraded.extend(forms); changed = True
                else:
                    upgraded.append(c)
        if changed:
            chosen = list(dict.fromkeys(upgraded))
            norm_method = "best_rhea_form"

        chosen = [c for c in dict.fromkeys(chosen) if c]
        # node id
        if len(chosen) == 1:
            node_id = f"CHEBI:{chosen[0].split(':')[-1]}"
        elif len(chosen) >= 2:
            node_id = _refmet_node(refmet, name)
        else:
            node_id = _refmet_node(refmet, name)
            if mapping_method != "unmapped":
                mapping_method = "unmapped"; confidence = "none"

        chebi_rhea = ";".join(chosen) if chosen else None
        in_human = any(c in ctx.human_rhea_chebi for c in chosen)
        n_rxn = sum(ctx.human_rxn_count.get(c, 0) for c in chosen)
        is_currency = any(c in ctx.currency for c in chosen)

        out_rows.append({
            "node_id": node_id,
            "metabolite": name,
            "refmet_id": refmet,
            "super_class": row.get("super_class"),
            "main_class": row.get("main_class"),
            "chebi_id": row.get("chebi_id"),
            "chebi_rhea": chebi_rhea,
            "normalization_method": norm_method,
            "mapping_method": mapping_method,
            "confidence": confidence,
            "is_complex_lipid": bool(is_lip),
            "is_currency": bool(is_currency),
            "in_human_rhea": bool(in_human),
            "n_human_rhea_reactions": int(n_rxn),
            "candidate_isomers": candidate_isomers,
        })

    out = pd.DataFrame(out_rows)

    # merge measured rows landing on the same CHEBI node -> keep both refmet_ids
    chebi_nodes = out[out.node_id.str.startswith("CHEBI:")]
    dup = chebi_nodes.groupby("node_id")["refmet_id"].nunique()
    for node in dup[dup > 1].index:
        ids = out[out.node_id == node]["refmet_id"].tolist()
        log.append(f"[node-merge] {node} <- refmet {ids}")

    # write unmapped to manual mappings (preserve existing)
    _append_manual(cfg, out, curation_dir, ctx)

    # outputs
    interim = _interim(cfg)
    out.to_parquet(os.path.join(interim, "metabolites_nonlipid.parquet"), index=False)
    _write_review(cfg, out, ctx)
    _write_report(cfg, out, ctx, obsolete_replacements, log)
    print(f"[mnet.metabolites] rows={len(out):,} mapped={ (out.mapping_method!='unmapped').sum():,} "
          f"unmapped={(out.mapping_method=='unmapped').sum():,}")
    return 0


def _append_manual(cfg, out, curation_dir, ctx):
    path = os.path.join(curation_dir, "manual_metabolite_mappings.csv")
    existing = pd.read_csv(path, dtype=str) if os.path.exists(path) else pd.DataFrame()
    have = set(existing["refmet_id"]) if "refmet_id" in existing.columns else set()
    unmapped = out[out.mapping_method == "unmapped"]
    new = []
    for _, r in unmapped.iterrows():
        if r["refmet_id"] in have:
            continue
        # best suggestion: exact name/inchikey attempt (already failed) -> leave blank
        new.append({"metabolite": r["metabolite"], "refmet_id": r["refmet_id"],
                    "suggested_chebi": "", "note": "auto: unmapped in Phase 2"})
    if new:
        cols = ["metabolite", "refmet_id", "suggested_chebi", "note"]
        merged = pd.concat([existing.reindex(columns=cols) if len(existing) else pd.DataFrame(columns=cols),
                            pd.DataFrame(new, columns=cols)], ignore_index=True)
        merged.to_csv(path, index=False)


def _write_review(cfg, out, ctx):
    rows = []
    for _, r in out.iterrows():
        final = r["chebi_rhea"]
        nm = "; ".join(ctx.chebi_name.get(c, "?") for c in str(final).split(";")) if final else ""
        rows.append({
            "metabolite": r["metabolite"], "refmet_id": r["refmet_id"],
            "original_chebi": r["chebi_id"], "final_chebi": final,
            "mapping_method": r["mapping_method"],
            "normalization_method": r["normalization_method"],
            "final_chebi_name": nm, "confidence": r["confidence"],
            "candidate_isomers": r["candidate_isomers"],
        })
    pd.DataFrame(rows).to_csv(
        os.path.join(cfg["paths"]["reports"], "mnet_phase2_mapping_review.csv"), index=False)


def _write_report(cfg, out, ctx, obsolete_replacements, log):
    from string_network.config import resolve_output_path  # noqa: F401
    mcfg = cfg.get("mnet", {}) or {}
    L = ["# mnet Phase 2 - metabolite mapping", ""]

    # Correction A: ge500 peek
    ppi_path = net.ppi_path(cfg)
    ppi = net.load_ppi(cfg)
    ppi_nodes = set(ppi["protein1"]) | set(ppi["protein2"])
    enz = _read(cfg, "rhea_enzymes")
    human = ctx._human_set()
    human_enz = set(enz[enz["uniprot"].isin(human)]["uniprot"])
    L.append("## Correction A - PPI network is the >= 500 build")
    L.append("")
    L.append(f"- ppi path: `{ppi_path}`")
    L.append(f"- rows: {len(ppi):,} | nodes: {len(ppi_nodes):,} | min score: {ppi['combined_score'].min():.1f}")
    L.append(f"- human Rhea enzymes: {len(human_enz):,}; of those already ge500 node ids: "
             f"**{len({u for u in human_enz if u in ppi_nodes}):,}**")
    L.append("")

    # Correction B: compound_type old vs new (from participants)
    parts = _read(cfg, "rhea_participants")
    L.append("## Correction B - compound_type from Rhea classes")
    L.append("")
    L.append("Rhea's compound class is now the source of truth (was ChEBI-formula heuristic).")
    L.append("")
    L.append("| compound_type | Phase 1 (formula, human uniq ChEBI) | Phase 2 (Rhea class, all uniq ChEBI) |")
    L.append("|---|---|---|")
    new_counts = parts.drop_duplicates("chebi_id").groupby("compound_type").size().to_dict()
    old_counts = {"small molecule": 3644, "generic": 624, "polymer": 2}
    for t in ["small molecule", "generic", "polymer", "other"]:
        L.append(f"| {t} | {old_counts.get(t,'-')} | {new_counts.get(t,0)} |")
    L.append(f"- participants with reactive-part/underlying ChEBI stored: "
             f"{parts['underlying_chebi'].notna().sum():,}")
    L.append("")

    # Correction C: measured currency / high-degree
    cmax = mcfg.get("currency_max_reactions", 150)
    meas_cur = out[out.is_currency]
    hi = out[out.n_human_rhea_reactions > cmax]
    L.append("## Correction C - measured currency / high-degree metabolites")
    L.append("")
    L.append(f"- measured metabolites flagged is_currency: {len(meas_cur)}")
    L.append(f"  {sorted(meas_cur['metabolite'].tolist())}")
    L.append(f"- measured metabolites in > {cmax} human reactions: "
             f"{sorted(set(hi['metabolite'].tolist()))}")
    L.append("- (no filtering applied yet - that is Phase 4)")
    L.append("")

    # coverage by super_class
    L.append("## Coverage by super_class")
    L.append("")
    L.append("| super_class | total | mapped | in human Rhea | unmapped |")
    L.append("|---|---|---|---|---|")
    g = out.copy()
    g["mapped"] = g.mapping_method != "unmapped"
    for sc, sub in g.groupby(g.super_class.fillna("(none)")):
        L.append(f"| {sc} | {len(sub)} | {int(sub.mapped.sum())} | "
                 f"{int(sub.in_human_rhea.sum())} | {int((~sub.mapped).sum())} |")
    L.append(f"| **TOTAL** | {len(g)} | {int(g.mapped.sum())} | "
             f"{int(g.in_human_rhea.sum())} | {int((~g.mapped).sum())} |")
    L.append("")

    L.append("## Counts per mapping_method / normalization_method")
    L.append("")
    L.append(f"- mapping_method: {out.mapping_method.value_counts().to_dict()}")
    L.append(f"- normalization_method: {out.normalization_method.value_counts().to_dict()}")
    L.append("")

    L.append("## Unmapped and low-confidence")
    L.append("")
    L.append(f"- unmapped ({(out.mapping_method=='unmapped').sum()}): "
             f"{sorted(out[out.mapping_method=='unmapped']['metabolite'].tolist())}")
    L.append(f"- low confidence ({(out.confidence=='low').sum()}): "
             f"{sorted(out[out.confidence=='low']['metabolite'].tolist())}")
    L.append("")

    L.append("## 10 charge-normalization examples (original -> Rhea)")
    L.append("")
    ex = out[(out.normalization_method == "ph7_3_mapping") & out.chebi_id.notna() &
             (out.chebi_id != out.chebi_rhea)].head(10)
    L.append("| metabolite | original | -> Rhea | Rhea name |")
    L.append("|---|---|---|---|")
    for _, r in ex.iterrows():
        L.append(f"| {r['metabolite']} | {r['chebi_id']} | {r['chebi_rhea']} | "
                 f"{ctx.chebi_name.get(r['chebi_rhea'],'?')} |")
    L.append("")
    L.append(f"- obsolete ChEBI replacements applied: {obsolete_replacements}")
    L.append(f"- log lines (ties / merges / obsolete): {len(log)} (see interim log)")
    L.append("")

    os.makedirs(cfg["paths"]["reports"], exist_ok=True)
    with open(os.path.join(cfg["paths"]["reports"], "mnet_phase2_metabolites.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
    with open(os.path.join(_interim(cfg), "mnet_phase2.log"), "w") as fh:
        fh.write("\n".join(log) + "\n")
    print(f"[mnet.metabolites] wrote reports/mnet_phase2_metabolites.md")
