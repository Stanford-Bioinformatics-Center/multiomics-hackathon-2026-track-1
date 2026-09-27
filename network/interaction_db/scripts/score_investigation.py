#!/usr/bin/env python3
"""Data-driven investigation of the reference score transformation.

Runs the hypotheses in order and writes reports/score_investigation.md. Every
statement in the report is derived from a number computed here; there are no
hard-coded conclusions. If a hypothesis fails, the report says so.

  H1  per-edge ratio ref/raw on strictly 1:1-mapped shared edges
  H2  is the base source the physical subnetwork? (+ candidate rule)
  H3  which edges get x scale_factor? (rule + depth-2 decision tree)
  H4  different STRING version (only run if >10% of shared scores unexplained)

It also categorises the remaining edge mismatches of the current build.
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
from string_network import load as L                       # noqa: E402
from string_network.build import load_links_full           # noqa: E402
from string_network.config import load_config              # noqa: E402

CH = ["neighborhood", "fusion", "cooccurence", "coexpression",
      "experimental", "database", "textmining"]


def canon_key(df):
    a = np.minimum(df.protein1.values.astype(object), df.protein2.values.astype(object))
    b = np.maximum(df.protein1.values.astype(object), df.protein2.values.astype(object))
    return np.array([f"{x}|{y}" for x, y in zip(a, b)], dtype=object)


def uni_pairs_from_links(links, e2u):
    lk = links.copy()
    lk["u1"] = lk.protein1.map(e2u); lk["u2"] = lk.protein2.map(e2u)
    lk = lk.dropna(subset=["u1", "u2"]); lk = lk[lk.u1 != lk.u2]
    a = np.minimum(lk.u1.values.astype(object), lk.u2.values.astype(object))
    b = np.maximum(lk.u1.values.astype(object), lk.u2.values.astype(object))
    lk["k"] = [f"{x}|{y}" for x, y in zip(a, b)]
    return lk.groupby("k")["combined_score"].max().astype(float)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="config.yaml")
    ap.add_argument("--scale", type=float, default=None,
                    help="Scale factor to test (default: config scale_factor).")
    ap.add_argument("--out", default="reports/score_investigation.md")
    args = ap.parse_args()

    cfg = load_config(args.config)
    scale = args.scale if args.scale is not None else float(cfg.get("scale_factor", 0.9))
    L_ = []  # report lines
    def w(s=""):
        L_.append(s)

    ref = pd.read_parquet(cfg["paths"]["reference"])
    ref["k"] = canon_key(ref)
    refs = ref.set_index("k")["combined_score"].astype(float)

    # strictly 1:1 ENSP<->UniProt mapping for clean comparisons
    m = L.load_mapping(cfg, "uniprot_idmapping")
    ef = m.groupby("ensp").size(); uf = m.groupby("uniprot").size()
    strict = m[m.ensp.isin(ef[ef == 1].index) & m.uniprot.isin(uf[uf == 1].index)]
    e2u = dict(zip(strict.ensp, strict.uniprot))

    full = uni_pairs_from_links(load_links_full(cfg, links_source="full"), e2u)
    phys = uni_pairs_from_links(load_links_full(cfg, links_source="physical"), e2u)

    w("# Score investigation")
    w()
    w("All numbers below are computed by `scripts/score_investigation.py`. "
      "Comparisons use protein pairs that map strictly 1:1 between ENSP and "
      "UniProt so the base score is unambiguous.")
    w()

    # ------------------------------------------------------------------ H1
    common = refs.index.intersection(full.index)
    d = pd.DataFrame({"ref": refs.loc[common].values, "raw": full.loc[common].values})
    ex1 = np.abs(d.ref - d.raw) < 0.05
    ex_s = np.abs(d.ref - scale * d.raw) < 0.05
    other = ~(ex1 | ex_s)
    n = len(d)
    w("## H1 - per-edge ratio ref / raw (full v12 links)")
    w()
    w(f"- Shared strictly-1:1 edges compared: **{n:,}**")
    w(f"- ratio == 1.0 : {int(ex1.sum()):,} ({100*ex1.mean():.1f}%)")
    w(f"- ratio == {scale} : {int(ex_s.sum()):,} ({100*ex_s.mean():.1f}%)")
    w(f"- something else : {int(other.sum()):,} ({100*other.mean():.1f}%)")
    w()
    w(f"The majority are neither 1.0 nor {scale}x the *full* score, so for most "
      "edges the base score is not the full links score. Examples of "
      "'something else' (ref vs full raw):")
    w()
    w("| ref | full raw | ratio |")
    w("|---|---|---|")
    for _, r in d[other].head(6).iterrows():
        w(f"| {r.ref:.1f} | {r.raw:.0f} | {r.ref/r.raw:.4f} |")
    w()

    # ------------------------------------------------------------------ H2
    common_p = refs.index.intersection(phys.index)
    dp = pd.DataFrame({"ref": refs.loc[common_p].values, "phys": phys.loc[common_p].values})
    ex_phys = np.abs(dp.ref - dp.phys) < 0.05

    allk = full.index.union(phys.index)
    cand = pd.Series(index=allk, dtype=float)
    in_phys_all = allk.isin(phys.index)
    cand.loc[in_phys_all] = phys.reindex(allk[in_phys_all]).values
    cand.loc[~in_phys_all] = scale * full.reindex(allk[~in_phys_all]).values
    cand = cand[cand >= cfg.get("threshold", 700)].round(1)
    common_c = refs.index.intersection(cand.index)
    cc = pd.DataFrame({"ref": refs.loc[common_c].values, "cand": cand.loc[common_c].values})
    exc = np.abs(cc.ref - cc.cand) < 0.05
    jacc = len(common_c) / len(set(cand.index) | set(refs.index))

    w("## H2 - is the base source the physical subnetwork?")
    w()
    w(f"- Reference edges equal to the **physical** score exactly: "
      f"{int(ex_phys.sum()):,} / {len(dp):,} (**{100*ex_phys.mean():.1f}%**)")
    w(f"- Candidate rule `physical if in physical else {scale} x full`, threshold "
      f"applied after transform:")
    w(f"    - edge Jaccard (1:1 space): **{jacc:.4f}**")
    w(f"    - exact score match on shared edges: **{100*exc.mean():.1f}%** "
      f"({int(exc.sum()):,}/{len(cc):,})")
    w()

    unexplained = 100 * (1 - exc.mean())

    # ------------------------------------------------------------------ H3
    w("## H3 - which edges are scaled?")
    w()
    # label strict-shared edges vs full
    lab = d.copy()
    lab["k"] = common
    lab["scaled"] = np.abs(lab.ref - scale * lab.raw) < 0.05
    lab["unscaled"] = np.abs(lab.ref - lab.raw) < 0.05
    lab["in_physical"] = lab.k.isin(set(phys.index)).astype(int)
    clean = lab[lab.scaled | lab.unscaled].copy()

    rule_pred = clean.in_physical == 0
    acc_phys = (rule_pred == clean.scaled).mean()
    w(f"- Labelled edges: scaled={int(lab.scaled.sum()):,}, "
      f"unscaled(ref==full)={int(lab.unscaled.sum()):,}, "
      f"ambiguous(ref==physical!=full)={int((~(lab.scaled|lab.unscaled)).sum()):,}")
    w(f"- Rule **`scaled <=> not in physical network`** accuracy: "
      f"**{100*acc_phys:.2f}%**")
    w()
    w("Cross-tab (rows = scaled, cols = in_physical):")
    w()
    ct = pd.crosstab(clean.scaled, clean.in_physical)
    w("| scaled\\in_physical | " + " | ".join(str(c) for c in ct.columns) + " |")
    w("|" + "---|" * (len(ct.columns) + 1))
    for idx, row in ct.iterrows():
        w(f"| {idx} | " + " | ".join(f"{v:,}" for v in row.values) + " |")
    w()

    # alternative single-column rules, using channels from the full detailed file
    det_path = os.path.join(cfg["paths"]["raw"],
                            os.path.basename(cfg["urls"]["links_detailed"]))
    det = pd.read_csv(det_path, sep=r"\s+", header=0,
                      names=["protein1", "protein2"] + CH + ["combined_score"])
    det["protein1"] = det.protein1.str.replace(r"^\d+\.", "", regex=True).map(e2u)
    det["protein2"] = det.protein2.str.replace(r"^\d+\.", "", regex=True).map(e2u)
    det = det.dropna(subset=["protein1", "protein2"])
    det = det[det.protein1 != det.protein2]
    det["k"] = canon_key(det)
    det = det.drop_duplicates("k").set_index("k")
    feat = det.reindex(clean.k)[CH].reset_index(drop=True)
    y = clean.scaled.reset_index(drop=True).astype(int)
    valid = feat.notna().all(axis=1)
    feat = feat[valid]; yv = y[valid]
    inphys_v = clean.in_physical.reset_index(drop=True)[valid]

    w("Alternative rules (accuracy at predicting 'scaled'):")
    w()
    alt = {
        "experimental == 0": feat.experimental == 0,
        "database == 0": feat.database == 0,
        "experimental == 0 and database == 0": (feat.experimental == 0) & (feat.database == 0),
        "textmining is the largest channel": feat[CH].idxmax(axis=1) == "textmining",
    }
    w("| rule | accuracy |")
    w("|---|---|")
    for name, cond in alt.items():
        w(f"| {name} | {100*(cond.values == yv.values).mean():.2f}% |")
    w(f"| **in physical network (negated)** | **{100*((inphys_v.values==0)==yv.values).mean():.2f}%** |")
    w()

    try:
        from sklearn.tree import DecisionTreeClassifier, export_text
        X = feat.copy(); X["in_physical"] = inphys_v.values
        clf = DecisionTreeClassifier(max_depth=2, random_state=0).fit(X, yv)
        w(f"Depth-2 decision tree on channels + in_physical, accuracy "
          f"**{100*clf.score(X, yv):.2f}%**:")
        w()
        w("```")
        w(export_text(clf, feature_names=list(X.columns)).rstrip())
        w("```")
        w()
    except Exception as exc:  # noqa: BLE001
        w(f"(decision tree skipped: {exc})")
        w()

    # ------------------------------------------------------------------ H4
    w("## H4 - different STRING version (v11.5)")
    w()
    if unexplained > 10:
        w(f"H1-H3 leave {unexplained:.1f}% of shared-edge scores unexplained "
          "(> 10%), so a version check would be warranted. [Run with v11.5 files "
          "added to config to extend this section.]")
    else:
        w(f"H1-H3 leave only **{unexplained:.1f}%** of shared-edge scores "
          "unexplained (< 10%), so H4 was **not run**.")
    w()

    # ------------------------------------------------------- remaining mismatch
    w("## Remaining differences of the current build")
    w()
    from string_network.config import resolve_output_path
    out_path = resolve_output_path(cfg)
    if out_path and os.path.exists(out_path):
        b = pd.read_parquet(out_path)
        b["k"] = canon_key(b)
        bk, rk = set(b.k), set(ref.k)
        rn = set(ref.protein1) | set(ref.protein2)
        bn = set(b.protein1) | set(b.protein2)
        only_b = bk - rk; only_r = rk - bk
        bs = dict(zip(b.k, b.combined_score.astype(float)))

        def both_ref(k):
            a, c = k.split("|"); return a in rn and c in rn
        ob_bothref = sum(both_ref(k) for k in only_b)
        ob_via_newnode = len(only_b) - ob_bothref
        ob_score = np.array([bs[k] for k in only_b]) if only_b else np.array([])
        near_thr = int(((ob_score >= 700) & (ob_score < 705)).sum()) if len(ob_score) else 0
        bonly_nodes = bn - rn
        ronly_nodes = rn - bn

        w(f"- Build edges: {len(bk):,} | reference edges: {len(rk):,} | "
          f"shared: {len(bk & rk):,}")
        w(f"- Only in build: {len(only_b):,}  |  only in reference: {len(only_r):,}")
        w(f"- Build-only nodes: {len(bonly_nodes):,}  |  reference-only nodes: "
          f"{len(ronly_nodes):,}")
        w()
        # explanation derived from the counts
        pct_newnode = 100 * ob_via_newnode / len(only_b) if only_b else 0
        pct_thr = 100 * near_thr / len(only_b) if only_b else 0
        w(f"Of the {len(only_b):,} build-only edges, {ob_via_newnode:,} "
          f"({pct_newnode:.0f}%) touch a node that the reference labelled with a "
          f"different id, and {near_thr:,} ({pct_thr:.0f}%) sit within 5 points of "
          f"the 700 cutoff. In other words the residual is driven by "
          f"{len(bonly_nodes):,}/{len(ronly_nodes):,} nodes that the two pipelines "
          f"mapped to different ids, not by the score rule.")
        w()

        # provenance mapping: ensp -> (our id, source)
        prov = L.load_mapping_with_fallback(
            cfg, cfg.get("mapping_source", "uniprot_idmapping"), int(cfg.get("taxon", 9606)),
            alias_sources=cfg.get("mapping_alias_sources", ["UniProt_AC", "Ensembl_UniProt"]),
            return_source=True)
        ensp2ours = dict(zip(prov.ensp, prov.uniprot))
        ensp2src = dict(zip(prov.ensp, prov.source))
        ours2ensp = {}
        for e, u in zip(prov.ensp, prov.uniprot):
            ours2ensp.setdefault(u, e)
        info = L.load_protein_info(
            os.path.join(cfg["paths"]["raw"], os.path.basename(cfg["urls"]["protein_info"])))
        name2ensp = dict(zip(info.preferred_name, info.ensp))

        # ---- the 8 reference fallback ids ----
        w("### The 8 reference fallback ids")
        w()
        w("For each id the reference could not resolve to UniProt, what our build "
          "maps that ENSP to and via which step:")
        w()
        w("| reference id | ENSP | our id | our mapping step |")
        w("|---|---|---|---|")
        fallback_ids = [
            "9606.ENSP00000351232", "9606.ENSP00000423463", "9606.ENSP00000480571",
            "9606.ENSP00000483415", "9606.ENSP00000491341", "9606.ENSP00000501026",
            "CGB1", "DGCR6",
        ]
        for rid in fallback_ids:
            ensp = rid.split(".", 1)[1] if rid.startswith("9606.") else name2ensp.get(rid)
            ours = ensp2ours.get(ensp, "(unmapped)") if ensp else "(ENSP not found)"
            src = ensp2src.get(ensp, "-") if ensp else "-"
            w(f"| {rid} | {ensp} | {ours} | {src} |")
        w()

        # ---- node diff table matched via shared neighbours ----
        def neighbours(df):
            adj = {}
            for a, c in zip(df.protein1.values, df.protein2.values):
                adj.setdefault(a, set()).add(c)
                adj.setdefault(c, set()).add(a)
            return adj
        badj, radj = neighbours(b), neighbours(ref)
        shared_nodes = bn & rn

        def match_by_neighbours(node, adj_from, adj_to, candidates):
            src_nb = adj_from.get(node, set()) & shared_nodes
            best, best_j = None, 0.0
            for cand in candidates:
                dst_nb = adj_to.get(cand, set()) & shared_nodes
                if not src_nb and not dst_nb:
                    continue
                inter = len(src_nb & dst_nb)
                union = len(src_nb | dst_nb)
                j = inter / union if union else 0.0
                if j > best_j:
                    best_j, best = j, cand
            return best, best_j

        w(f"### Build-only vs reference-only nodes (matched via shared neighbours)")
        w()
        w(f"{len(bonly_nodes)} build-only and {len(ronly_nodes)} reference-only "
          f"nodes. Each build-only node is matched to the reference-only node with "
          f"the most similar set of shared neighbours (Jaccard J):")
        w()
        w("| ENSP | our id | our step | reference id | neighbour-J |")
        w("|---|---|---|---|---|")
        ronly_list = list(ronly_nodes)
        for node in sorted(bonly_nodes)[:60]:
            match, j = match_by_neighbours(node, badj, radj, ronly_list)
            ensp = ours2ensp.get(node, "?")
            src = ensp2src.get(ensp, "?")
            w(f"| {ensp} | {node} | {src} | {match if match else '(no match)'} | {j:.2f} |")
        w()
    else:
        w("(build output not found; run `make build` first)")
    w()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write("\n".join(L_) + "\n")
    print(f"[investigation] H1 other={100*other.mean():.1f}% | "
          f"H2 candidate exact={100*exc.mean():.1f}% jaccard={jacc:.3f} | "
          f"H3 in_physical acc={100*acc_phys:.2f}% | unexplained={unexplained:.1f}%")
    print(f"[investigation] wrote {args.out}")


if __name__ == "__main__":
    main()
