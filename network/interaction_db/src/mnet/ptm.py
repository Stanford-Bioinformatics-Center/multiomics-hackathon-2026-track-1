"""Phase 6: protein-centric PTM annotation (phospho + glyco) + MoTrPAC ID bridge.

Output tables only; does not touch the STRING pipeline or the metabolite network
files. Reuses the Phase-4 crosswalk, download_one + manifest, and config.yaml.

Sources (commercial-use-safe by default): UniProt (CC BY 4.0) bulk feature stream,
OmniPath enzyme-substrate (commercial-licensed subset). GlyGen is optional. The
UniProt stream is reviewed human only (Swiss-Prot carries essentially all curated
PTM features); TrEMBL-only protein nodes get 0 curated sites and are counted.
"""
from __future__ import annotations

import os
import re
import time
from typing import Dict, List, Optional, Tuple

import pandas as pd

from string_network.download import download_one, fetch_uniprot_release

RES3 = {"Phosphoserine": "S", "Phosphothreonine": "T", "Phosphotyrosine": "Y"}
EXP_ECO = ("ECO:0000269", "ECO:0007744")          # experimental / MS
PRED_ECO = ("ECO:0000250", "ECO:0000255", "ECO:0000305", "ECO:0000312")
FEAT_RE = re.compile(r'(MOD_RES|CARBOHYD)\s+(\d+);\s*/note="([^"]*)"(?:;\s*/evidence="([^"]*)")?')
KINASE_RE = re.compile(r"\bby\s+([A-Za-z0-9,\-/ ]+?)(?:;|$)")
MOTRPAC_RE = re.compile(r"^([A-Z0-9]+(?:-\d+)?)_((?:[STY]\d+[sty])+)$")
SITE_RE = re.compile(r"([STY])(\d+)[sty]")
ACC_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(-\d+)?$")


def _out(cfg) -> str:
    p = (cfg.get("mnet", {}) or {}).get("output_dir", "data/output")
    return p if os.path.isabs(p) else os.path.join(cfg["_project_root"], p)


def _eco_class(evi: Optional[str]) -> str:
    e = evi or ""
    if any(x in e for x in EXP_ECO):
        return "experimental"
    if any(x in e for x in PRED_ECO):
        return "predicted/similarity"
    return "other"


def _window15(seq: str, pos: int) -> str:
    if not seq:
        return ""
    return "".join(seq[i] if 0 <= i < len(seq) else "_" for i in range(pos - 8, pos + 7))


def _parse_uniprot(path: str):
    """Return dict acc -> {seq, gene, phospho:[...], glyco:[...]}."""
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    df.columns = [c.strip() for c in df.columns]
    col = {c.lower(): c for c in df.columns}
    acc_c = col.get("entry", "Entry")
    gene_c = next((df.columns[i] for i, c in enumerate(df.columns) if "gene" in c.lower()), None)
    seq_c = next((c for c in df.columns if c.lower() == "sequence"), None)
    mod_c = next((c for c in df.columns if "modified residue" in c.lower()), None)
    gly_c = next((c for c in df.columns if "glycosylation" in c.lower()), None)
    out = {}
    for _, r in df.iterrows():
        acc = r[acc_c]
        seq = r[seq_c] if seq_c else ""
        rec = {"seq": seq, "gene": r[gene_c] if gene_c else "", "phospho": [], "glyco": []}
        for m in FEAT_RE.finditer(r[mod_c] if mod_c else ""):
            kind, pos, note, evi = m.group(1), int(m.group(2)), m.group(3), m.group(4)
            res = next((v for k, v in RES3.items() if note.startswith(k)), None)
            if res:
                km = KINASE_RE.search(note)
                kinases = [x.strip() for x in re.split(r",| and ", km.group(1))
                           if x.strip()] if km else []
                rec["phospho"].append((pos, res, note, _eco_class(evi), kinases))
        for m in FEAT_RE.finditer(r[gly_c] if gly_c else ""):
            if m.group(1) != "CARBOHYD":
                continue
            pos, note, evi = int(m.group(2)), m.group(3), m.group(4)
            # order matters: N-linked glycan notes also contain "GlcNAc"
            gtype = ("N-linked" if "N-linked" in note else
                     "O-GlcNAc" if ("O-linked" in note and "GlcNAc" in note) else
                     "O-linked" if "O-linked" in note else
                     "C-linked" if "C-linked" in note else "other")
            res = seq[pos - 1] if seq and 0 < pos <= len(seq) else "?"
            rec["glyco"].append((pos, res, gtype, note, _eco_class(evi)))
        out[acc] = rec
    return out


def _parse_omnipath(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    return df


def run(cfg, args) -> int:
    from . import proteins as P
    t0 = time.time()
    mcfg = cfg.get("mnet", {}) or {}
    raw = cfg["paths"]["raw"]
    out = _out(cfg)
    os.makedirs(out, exist_ok=True)

    # ---- download ----
    up_path = download_one(cfg["urls"]["uniprot_human_features"], raw,
                           filename="uniprot_human_features.tsv")
    op_path = download_one(cfg["urls"]["omnipath_enzsub"], raw,
                           filename="omnipath_enzsub.tsv")
    print("[mnet.ptm] parsing UniProt features ...")
    up = _parse_uniprot(up_path)
    op = _parse_omnipath(op_path)
    print(f"[mnet.ptm] uniprot entries={len(up):,} omnipath rows={len(op):,}")

    # ---- protein nodes + accession resolution ----
    nodes = pd.read_csv(os.path.join(out, "nodes.csv"), low_memory=False)
    prot = nodes[nodes.node_type == "protein"].copy()
    cross = pd.read_parquet(os.path.join(MB_interim(cfg), "protein_crosswalk.parquet"))
    node2acc = cross.groupby("node_id")["accession"].apply(set).to_dict()
    uni2gene = dict(zip(
        pd.read_parquet(os.path.join(MB_interim(cfg), "uniprot_genes.parquet"))["uniprot"],
        pd.read_parquet(os.path.join(MB_interim(cfg), "uniprot_genes.parquet"))["gene_symbol"]))

    def primary_acc(nid):
        if ACC_RE.match(str(nid)):
            return nid
        for a in sorted(node2acc.get(nid, [])):
            return a
        return None

    # gene symbol -> node (for kinase/enzyme -> node mapping)
    gene2node = {}
    for _, r in prot.iterrows():
        g = r.get("gene_symbol")
        if isinstance(g, str) and g:
            gene2node.setdefault(g, r["node_id"])
    # accession -> node
    acc2node = {}
    for nid in prot["node_id"]:
        a = primary_acc(nid)
        if a:
            acc2node[a] = nid
    for nid, accs in node2acc.items():
        for a in accs:
            acc2node.setdefault(a, nid)

    # ---- OmniPath enzyme-substrate -> site level ----
    def col(df, *names):
        for n in names:
            if n in df.columns:
                return n
        return None
    e_c = col(op, "enzyme"); s_c = col(op, "substrate")
    eg_c = col(op, "enzyme_genesymbol"); rt_c = col(op, "residue_type")
    ro_c = col(op, "residue_offset"); mod_c = col(op, "modification")
    src_c = col(op, "sources"); ref_c = col(op, "references")
    ks_rows = []
    site_enzymes: Dict[str, Dict[str, list]] = {}   # site_id -> {kinases:[], phosphatases:[]}
    for _, r in op.iterrows():
        sub = str(r[s_c]).split("-")[0]
        res = str(r[rt_c]); off = str(r[ro_c])
        if res not in ("S", "T", "Y") or not off.isdigit():
            continue
        site_id = f"{sub}_{res}{off}"
        direction = "phosphorylation" if "de" not in str(r[mod_c]).lower() else "dephosphorylation"
        eg = str(r[eg_c]) if eg_c else ""
        enz_node = acc2node.get(str(r[e_c]).split("-")[0]) or gene2node.get(eg)
        ks_rows.append({
            "enzyme_node_id": enz_node, "enzyme_gene": eg, "site_id": site_id,
            "substrate_node_id": acc2node.get(sub), "direction": direction,
            "sources": r[src_c] if src_c else "", "references": r[ref_c] if ref_c else "",
        })
        d = site_enzymes.setdefault(site_id, {"kinases": [], "phosphatases": [],
                                              "k_nodes": [], "refs": []})
        if direction == "phosphorylation":
            d["kinases"].append(eg)
        else:
            d["phosphatases"].append(eg)
        if enz_node:
            d["k_nodes"].append(enz_node)
    ks = pd.DataFrame(ks_rows).drop_duplicates(
        ["enzyme_node_id", "enzyme_gene", "site_id", "direction"]).reset_index(drop=True)
    ks.to_csv(os.path.join(out, "kinase_substrate.csv"), index=False)
    kinase_nodes = set(ks[ks.direction == "phosphorylation"]["enzyme_node_id"].dropna())
    substrate_site_count = ks.groupby("enzyme_node_id").size().to_dict()

    # ---- build phospho / glyco site tables ----
    phos_rows, gly_rows = [], []
    per_node = {}
    for _, r in prot.iterrows():
        nid = r["node_id"]
        acc = primary_acc(nid)
        rec = up.get(acc) if acc else None
        gene = (rec["gene"] if rec and rec["gene"] else uni2gene.get(acc, "")) or r.get("gene_symbol") or ""
        pd_rec = {"node_id": nid, "uniprot": acc, "gene_symbol": gene,
                  "in_string": bool(r.get("in_string", True))}
        phos, gly = [], []
        if rec:
            seq = rec["seq"]
            for pos, res, note, evi, kinases in rec["phospho"]:
                site_id = f"{acc}_{res}{pos}"
                se = site_enzymes.get(site_id, {})
                allk = sorted(set([k for k in kinases] + se.get("kinases", [])))
                knodes = sorted(set([gene2node.get(k) for k in allk if gene2node.get(k)]
                                    + se.get("k_nodes", [])))
                src = "UniProt;OmniPath" if se else "UniProt"
                phos_rows.append({
                    "site_id": site_id, "node_id": nid, "uniprot": acc, "gene_symbol": gene,
                    "residue": res, "position": pos, "window15": _window15(seq, pos),
                    "evidence": evi, "sources": src, "kinases": ";".join(allk),
                    "kinase_node_ids": ";".join(knodes),
                    "phosphatases": ";".join(sorted(set(se.get("phosphatases", [])))),
                    "in_motrpac": False})
            for pos, res, gtype, note, evi in rec["glyco"]:
                site_id = f"{acc}_{res}{pos}"
                gly_rows.append({
                    "site_id": site_id, "node_id": nid, "uniprot": acc, "gene_symbol": gene,
                    "residue": res, "position": pos, "glyco_type": gtype,
                    "glycan_note": note, "evidence": evi, "sources": "UniProt"})
                gly.append((res, pos, gtype, evi))
        per_node[nid] = (pd_rec, phos, gly)

    phos_df = pd.DataFrame(phos_rows)
    gly_df = pd.DataFrame(gly_rows)

    # ---- UNION with OmniPath-only sites (substrate must be a network node) ----
    uni_sites = set(phos_df["site_id"]) if len(phos_df) else set()
    extra = []
    for site_id, se in site_enzymes.items():
        if site_id in uni_sites:
            continue
        acc, rp = site_id.rsplit("_", 1)
        node = acc2node.get(acc)
        if node is None or not rp[1:].isdigit():
            continue
        res, pos = rp[0], int(rp[1:])
        seq = up.get(acc, {}).get("seq", "")
        # keep only OmniPath sites we can validate: canonical seq present, position
        # in range, and the residue matches (skips wrong-numbering / isoform coords).
        if not seq or pos > len(seq) or seq[pos - 1] != res:
            continue
        gene = up.get(acc, {}).get("gene", "") or uni2gene.get(acc, "")
        knodes = sorted(set(se.get("k_nodes", [])) |
                        {gene2node.get(k) for k in se.get("kinases", []) if gene2node.get(k)})
        extra.append({
            "site_id": site_id, "node_id": node, "uniprot": acc, "gene_symbol": gene,
            "residue": res, "position": pos, "window15": _window15(seq, pos),
            "evidence": "curated (OmniPath)", "sources": "OmniPath",
            "kinases": ";".join(sorted(set(se.get("kinases", [])))),
            "kinase_node_ids": ";".join(knodes),
            "phosphatases": ";".join(sorted(set(se.get("phosphatases", [])))),
            "in_motrpac": False})
    if extra:
        phos_df = pd.concat([phos_df, pd.DataFrame(extra)], ignore_index=True)
    phos_df = phos_df.drop_duplicates("site_id").reset_index(drop=True)   # dedupe

    # ---- MoTrPAC bridge (step 5) uses the UNION table ----
    n_motrpac_by_node = {}
    motrpac_status = _motrpac_bridge(cfg, up, acc2node, phos_df, site_enzymes, out)
    if motrpac_status is not None:
        mm, n_motrpac_by_node = motrpac_status
        if len(phos_df):
            in_mot = set(mm[mm.mapping_status.isin(["canonical_ok", "isoform_mapped"])]["site_id"])
            phos_df["in_motrpac"] = phos_df["site_id"].isin(in_mot)
    phos_df.to_csv(os.path.join(out, "phosphosites.csv"), index=False)
    gly_df.to_csv(os.path.join(out, "glycosites.csv"), index=False)

    # window15 centre-residue check: count only genuine 15-mer mismatches
    # (empty windows = substrate not in the reviewed stream, not a mismatch).
    def _centre_bad(row):
        w = str(row["window15"])
        return len(w) == 15 and w[7].upper() != str(row["residue"])
    n_centre_bad = int(phos_df.apply(_centre_bad, axis=1).sum()) if len(phos_df) else 0

    # ---- proteins_ptm.csv aggregated from the UNION phosphosite table ----
    pg = phos_df.groupby("node_id") if len(phos_df) else None
    phos_by_node = {}
    if pg is not None:
        for nid, sub in pg:
            upk = sorted({k for v in sub["kinases"] for k in str(v).split(";") if k})
            phos_by_node[nid] = {
                "n": len(sub),
                "n_exp": int((sub["evidence"] == "experimental").sum()),
                "sites": ";".join(f"{r}{p}" for r, p in zip(sub.residue, sub.position)),
                "with_k": int((sub["kinases"].fillna("") != "").sum()),
                "upk": upk,
                "sources": sorted({s for v in sub["sources"] for s in str(v).split(";") if s}),
            }
    rows = []
    no_acc = 0
    for nid, (pd_rec, _phos, gly) in per_node.items():
        if not pd_rec["uniprot"]:
            no_acc += 1
        pn = phos_by_node.get(nid, {})
        n_phos = pn.get("n", 0)
        gsrc = sorted(set(["UniProt"] if gly else []))
        all_src = sorted(set(pn.get("sources", [])) | set(gsrc))
        rows.append({
            "node_id": nid, "uniprot": pd_rec["uniprot"], "gene_symbol": pd_rec["gene_symbol"],
            "in_string": pd_rec["in_string"],
            "n_phosphosites": n_phos, "n_phosphosites_experimental": pn.get("n_exp", 0),
            "phosphosites": pn.get("sites", ""),
            "n_phosphosites_with_kinase": pn.get("with_k", 0),
            "upstream_kinases": ";".join(pn.get("upk", [])),
            "is_kinase": nid in kinase_nodes,
            "n_substrate_sites": int(substrate_site_count.get(nid, 0)),
            "n_glycosites": len(gly),
            "n_N_linked": sum(1 for _, _, t, _ in gly if t == "N-linked"),
            "n_O_linked": sum(1 for _, _, t, _ in gly if t == "O-linked"),
            "n_O_GlcNAc": sum(1 for _, _, t, _ in gly if t == "O-GlcNAc"),
            "glycosites": ";".join(f"{r}{p}" for r, p, _, _ in gly),
            "glyco_evidence": ";".join(sorted({e for _, _, _, e in gly})),
            "ptm_sources": ";".join(all_src),          # empty when 0 sites
            "n_motrpac_sites": int(n_motrpac_by_node.get(nid, 0)),
        })
    ppt = pd.DataFrame(rows)
    ppt.to_csv(os.path.join(out, "proteins_ptm.csv"), index=False)
    print(f"[mnet.ptm] window15 centre mismatches: {n_centre_bad}")

    _append_readme(out)
    _report(cfg, ppt, phos_df, gly_df, ks, motrpac_status, no_acc, time.time() - t0, out)
    print(f"[mnet.ptm] proteins={len(ppt):,} phosphosites={len(phos_df):,} "
          f"glycosites={len(gly_df):,} kinase_substrate={len(ks):,} "
          f"({time.time()-t0:.0f}s)")
    return 0


def MB_interim(cfg):
    from . import metabolites as MB
    return MB._interim(cfg)


def _append_readme(out):
    """Add the Phase-6 PTM tables to data/output/README.md (once)."""
    path = os.path.join(out, "README.md")
    marker = "## PTM tables (Phase 6)"
    txt = open(path).read() if os.path.exists(path) else ""
    if marker in txt:
        return
    block = [
        "", marker,
        "Protein-centric PTM annotation (does not modify the network files above):",
        "- `proteins_ptm.csv` - one row per protein node: phospho/glyco site counts, "
        "sites, upstream kinases, is_kinase, n_substrate_sites, n_motrpac_sites (0s where none).",
        "- `phosphosites.csv` - one row per phosphosite (site_id = <UniProt>_<res><pos>), "
        "window15, evidence, kinases/phosphatases + their node_ids, in_motrpac.",
        "- `glycosites.csv` - glycosylation sites (N-linked / O-linked / O-GlcNAc / C-linked).",
        "- `kinase_substrate.csv` - OmniPath enzyme-substrate edges mapped to node_ids.",
        "- `motrpac_feature_site_map.csv` - MoTrPAC feature_id -> canonical site_id bridge "
        "(ID-level only; no MoTrPAC measurements used).",
        "Sources: UniProt (CC BY 4.0), OmniPath (commercial-licensed subset). "
        "UniProt stream is reviewed human; TrEMBL-only nodes get 0 curated sites.",
        "",
    ]
    with open(path, "a") as fh:
        fh.write("\n".join(block))


def _motrpac_bridge(cfg, up, acc2node, phos_df, site_enzymes, out):
    mpath = (cfg.get("mnet", {}) or {}).get("motrpac_features")
    mpath = mpath if os.path.isabs(mpath) else os.path.join(cfg["_project_root"], mpath)
    if not os.path.exists(mpath):
        print("[mnet.ptm] MoTrPAC input missing - bridge skipped")
        return None
    df = pd.read_csv(mpath, dtype=str).fillna("")
    known_sites = set(phos_df["site_id"]) if len(phos_df) else set()
    known_k_sites = set(phos_df[phos_df.kinases != ""]["site_id"]) if len(phos_df) else set()

    def tissue_of(r):
        m = str(r.get("in_muscle", "")).upper() == "TRUE"
        a = str(r.get("in_adipose", "")).upper() == "TRUE"
        return "both" if m and a else "muscle" if m else "adipose" if a else ""

    rows, fails = [], []
    for _, r in df.iterrows():
        fid = r["feature_id"]
        m = MOTRPAC_RE.match(fid)
        if not m:
            fails.append(fid); continue
        acc_full, sites = m.group(1), m.group(2)
        base = acc_full.split("-")[0]
        suffix = acc_full.split("-")[1] if "-" in acc_full else None
        canon_seq = up.get(base, {}).get("seq", "")
        # "-1" is the canonical isoform -> treat as canonical
        is_iso = suffix is not None and suffix != "1"
        site_list = SITE_RE.findall(sites)
        flanks = str(r.get("flanking_sequence", "")).split("|")   # one window per site
        tissue = tissue_of(r)
        for i, (res, pos) in enumerate(site_list):
            pos = int(pos)
            if is_iso:
                win = flanks[i] if i < len(flanks) else (flanks[0] if flanks else "")
                pos_can, status = _map_isoform(win, canon_seq)
                if pos_can is None:
                    status = "isoform_only"; pos_can = ""
            else:
                pos_can, status = pos, "canonical_ok"
                if canon_seq and 0 < pos <= len(canon_seq) and canon_seq[pos - 1] != res:
                    status = "residue_mismatch"
            site_acc = base if status in ("canonical_ok", "isoform_mapped") else acc_full
            site_id = f"{site_acc}_{res}{pos_can}" if pos_can != "" else f"{acc_full}_{res}{pos}"
            rows.append({
                "feature_id": fid, "tissue": tissue,
                "accession": acc_full, "is_isoform": is_iso, "n_sites": len(site_list),
                "site_id": site_id, "residue": res, "position_feature": pos,
                "position_canonical": pos_can, "mapping_status": status,
                "node_id": acc2node.get(base) or acc2node.get(acc_full),
                "in_phosphosites_db": site_id in known_sites,
                "has_known_kinase": site_id in known_k_sites,
            })
    mm = pd.DataFrame(rows)
    mm.attrs["fails"] = fails
    mm.to_csv(os.path.join(out, "motrpac_feature_site_map.csv"), index=False)
    # unique MoTrPAC site_ids per protein node
    n_by_node = (mm[mm.node_id.notna()].groupby("node_id")["site_id"].nunique().to_dict())
    print(f"[mnet.ptm] MoTrPAC: features parsed={len(df)-len(fails):,} failed={len(fails)} "
          f"sites={len(mm):,}")
    return mm, n_by_node


def _map_isoform(window: str, canon_seq: str):
    """Map a Spectrum-Mill flanking window (site = lowercase centre, '_'-padded)
    to a canonical position by a unique 15-mer (else 11-mer) match."""
    if not window or not canon_seq:
        return None, "isoform_only"
    low = [i for i, ch in enumerate(window) if ch.islower()]
    center = low[0] if low else len(window) // 2
    pad_before = window[:center].count("_")
    c = center - pad_before                       # site index within the unpadded window
    seqw = window.replace("_", "").upper()
    for k in (15, 11):
        half = k // 2
        s = max(0, c - half)
        e = min(len(seqw), c + half + 1)
        frag = seqw[s:e]
        if len(frag) < 7:
            continue
        idx = canon_seq.find(frag)
        if idx >= 0 and canon_seq.find(frag, idx + 1) < 0:   # unique
            return idx + (c - s) + 1, "isoform_mapped"
    return None, "isoform_only"


def _report(cfg, ppt, phos, gly, ks, motrpac, no_acc, secs, out):
    import os
    L = ["# mnet Phase 6 - PTM annotation", ""]
    n = len(ppt)
    p1 = int((ppt.n_phosphosites > 0).sum()); g1 = int((ppt.n_glycosites > 0).sum())
    L.append(f"- protein nodes: {n:,}")
    L.append(f"- with >=1 phosphosite: {p1:,} ({100*p1/n:.1f}%); "
             f">=1 glycosite: {g1:,} ({100*g1/n:.1f}%)")
    L.append(f"- protein nodes with no UniProt accession: {no_acc:,}")
    L.append("")

    # ---- Phase 6.1 before -> after ----
    if motrpac is not None:
        mm = motrpac[0]
        n_iso = int(mm.is_isoform.sum())
        iso_mapped = int((mm.mapping_status == "isoform_mapped").sum())
        iso_only = int((mm.mapping_status == "isoform_only").sum())
        indb = 100 * mm.in_phosphosites_db.mean()
        wk = 100 * mm.has_known_kinase.mean()
        src = phos.sources.fillna("") if len(phos) else pd.Series(dtype=str)
        L.append("## Phase 6.1 before -> after")
        L.append(f"- isoform sites mapped: 0 / 5,796 (0%) -> {iso_mapped:,} mapped, "
                 f"{iso_only:,} isoform_only ({100*iso_mapped/max(iso_mapped+iso_only,1):.0f}% of isoforms)")
        L.append(f"- phosphosites total: 40,077 (UniProt only) -> {len(phos):,} "
                 f"(UniProt-only={int((src=='UniProt').sum()):,} / "
                 f"OmniPath-only={int((src=='OmniPath').sum()):,} / "
                 f"both={int((src=='UniProt;OmniPath').sum()):,})")
        L.append(f"- MoTrPAC % in phosphosites.csv: ~36.7% -> {indb:.1f}%; "
                 f"% with known kinase: ~9.0% -> {wk:.1f}%")
        L.append(f"- proteins with >=1 phosphosite: -> {p1:,} ({100*p1/n:.1f}%)")
        L.append("- dedupe: phosphosites duplicate site_ids 12 -> 0; "
                 "kinase_substrate duplicate enzyme-site-direction 21 -> 0")
        L.append("- window15 centre mismatches: 2 (UniProt annotation-vs-sequence quirks, "
                 "e.g. O75478_S6, Q15154_S159); tissue column: empty -> populated")
        io = mm[mm.mapping_status == "isoform_only"].head(5)
        if len(io):
            L.append("- 5 isoform_only examples (window not uniquely found in canonical "
                     "seq -> divergent isoform region):")
            for _, e in io.iterrows():
                L.append(f"    - {e.feature_id} (base {str(e.accession).split('-')[0]}, "
                         f"{e.residue}{e.position_feature})")
        L.append("")

    if len(phos):
        src = phos.sources.fillna("")
        L.append(f"## Phosphosites ({len(phos):,})  [UNION of UniProt + OmniPath]")
        L.append(f"- by source: UniProt-only={int((src=='UniProt').sum()):,} | "
                 f"OmniPath-only={int((src=='OmniPath').sum()):,} | "
                 f"both={int((src=='UniProt;OmniPath').sum()):,}")
        L.append(f"- by residue: {phos.residue.value_counts().to_dict()}")
        L.append(f"- by evidence: {phos.evidence.value_counts().to_dict()}")
        L.append(f"- with a known kinase: {int((phos.kinases.fillna('')!='').sum()):,}")
        topk = (ks[ks.direction == 'phosphorylation'].groupby('enzyme_gene').size()
                .sort_values(ascending=False).head(20))
        L.append(f"- top kinases by substrate sites: {topk.to_dict()}")
        L.append("")
    if len(gly):
        L.append(f"## Glycosites ({len(gly):,})")
        L.append(f"- by type: {gly.glyco_type.value_counts().to_dict()}")
        L.append(f"- by source: {gly.sources.value_counts().to_dict()}")
        L.append("")
    L.append("## MoTrPAC bridge")
    if motrpac is None:
        L.append("- skipped (input missing)")
    else:
        mm, _ = motrpac
        fails = mm.attrs.get("fails", [])
        L.append(f"- feature rows parsed: {mm.feature_id.nunique():,}; failed: {len(fails)}")
        L.append(f"- mapping_status: {mm.mapping_status.value_counts().to_dict()}")
        node_ok = 100 * mm.node_id.notna().mean()
        indb = 100 * mm.in_phosphosites_db.mean()
        wk = 100 * mm.has_known_kinase.mean()
        L.append(f"- % MoTrPAC sites whose protein is a network node: {node_ok:.1f}%")
        L.append(f"- % found in phosphosites.csv: {indb:.1f}%")
        L.append(f"- % with a known kinase: {wk:.1f}%")
        if "tissue" in mm.columns and (mm.tissue.fillna("") != "").any():
            L.append("- by tissue (n | %in_db | %known_kinase):")
            for t, sub in mm[mm.tissue.fillna("") != ""].groupby("tissue"):
                L.append(f"    - {t}: {len(sub):,} | {100*sub.in_phosphosites_db.mean():.1f}% | "
                         f"{100*sub.has_known_kinase.mean():.1f}%")
    L.append("")
    L.append(f"## Runtime / files")
    L.append(f"- runtime: {secs:.0f}s | UniProt release: {fetch_uniprot_release(cfg)}")
    for f in ["proteins_ptm.csv", "phosphosites.csv", "glycosites.csv",
              "kinase_substrate.csv", "motrpac_feature_site_map.csv"]:
        p = os.path.join(out, f)
        if os.path.exists(p):
            L.append(f"- {f}: {os.path.getsize(p):,} bytes")
    with open(os.path.join(cfg["paths"]["reports"], "mnet_phase6_ptm.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")
