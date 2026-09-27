"""The start-node menu: the nodes on the team's three story slides that make good walk starts.

The three stories are step 20's primary results (20_stories.csv) and the team's story slides
(story_figures_FINAL_2026-09-27/clean_for_slides: muscle_t2d.png, blood_t2d.png, blood_ageing.png):
  T2D muscle    Öhman 2021 T2D muscle proteins; ENDURANCE reverses them (the endurance arm's muscle response)
  T2D blood     UK Biobank incident-T2D plasma proteins (Gadd 2024); RESISTANCE reverses them (blood protein response)
  ageing blood  UK Biobank ageing plasma proteins (Sun 2023); RESISTANCE reverses them (blood protein response)
The menu offers ONLY nodes labelled on those slides (exvideo/slide_nodes.csv, transcribed from the three PNGs on
2026-09-27; the slides' own branch, story-networks @ 50886d1, is not pushed) that are ROBUST starts: figure 17 nodes from
which every 3-step walk (never revisiting a node) goes the full distance, whichever way it turns (walk.robust_starts).
For each, the story's disease z (20_disease_scores.csv, if the node is in the set) and the winning arm's response
(mean of the story tissue's cells in 01_nodes_EE / RE.csv), exactly as step 20 computes them; the reversal score
-(z / sd z) x (response / sd response) over the set's altered proteins orders the menu (> 0 = the arm pushes it back).
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path
from typing import Dict, List, Optional, Sequence

from .errors import InputError
from .network import Network
from .walk import robust_starts

SLIDE_FILE = Path(__file__).resolve().parent / "slide_nodes.csv"

STORIES = [("T2D muscle", "ohman_2021", "muscle", "EE"),
           ("T2D blood", "gadd_2024_ukb_incident_T2D", "blood_prot", "RE"),
           ("ageing blood", "sun_2023_ukb_age", "blood_prot", "RE")]
ARM_NAME = {"EE": "endurance", "RE": "resistance"}


def _rows(path: Path) -> List[Dict[str, str]]:
    if not path.exists():
        raise InputError(f"missing pipeline output {path} (run: bash network/run_all.sh; steps 20 and 01 make it)")
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _tissue_mean(rows: List[Dict[str, str]], tissue: str) -> Dict[str, float]:
    out = {}
    for r in rows:
        v = [float(r[k]) for k in r if k.startswith(tissue + "_") and r[k] not in ("", "NA")]
        if v:
            out[r["gene_symbol"]] = sum(v) / len(v)
    return out


def slide_nodes() -> Dict[str, List[str]]:
    """The nodes labelled on each story slide (exvideo/slide_nodes.csv), by story."""
    out: Dict[str, List[str]] = {}
    for r in _rows(SLIDE_FILE):
        out.setdefault(r["story"], []).append(r["node"])
    return out


def scope_nodes(story: str, net: Network) -> List[str]:
    """The figure 17 nodes on one story slide: the walk's whole world when the walk is scoped to that slide."""
    return sorted(n for n in slide_nodes().get(story, []) if n in net.nodes)


def story_markers(out: Path, net: Network, per_story: int = 99, stories: Optional[Sequence[str]] = None) -> List[Dict[str, object]]:
    """Each story slide's robust-start nodes, best reversal first (nodes outside the story's altered set last). Robust
    = every 3-step walk that stays ON THAT SLIDE goes the full distance (the walk is confined to the slide)."""
    from .network import subnetwork
    slides = slide_nodes()
    meta = {r["set"]: r for r in _rows(out / "20_disease_sets.csv")}
    ds = _rows(out / "20_disease_scores.csv")
    resp_rows = {a: _rows(out / f"01_nodes_{a}.csv") for a in ("EE", "RE")}
    # red/blue links: the arm-specific edges the video draws bold (the figure 17 rule, as in network.build_facts)
    pool = sorted(abs(w) for e in net.edges.values() for w in (e.w_ee, e.w_re)); tau = pool[int(0.75 * (len(pool) - 1))] if pool else 0.0
    red_blue = lambda n: sum(1 for e in net.neighbours(n) if (e.w_ee >= tau) != (e.w_re >= tau))
    picked = []
    for story, s, tissue, arm in STORIES:
        if stories is not None and story not in stories:
            continue
        ok = set(robust_starts(subnetwork(net, slides.get(story, []))))
        if s not in meta:
            raise InputError(f"{s} is not in 20_disease_sets.csv (re-run step 20)")
        alpha, sig_only = float(meta[s]["alpha"]), meta[s]["sig_only"] == "TRUE"
        resp = _tissue_mean(resp_rows[arm], tissue)
        zs = {r["gene"]: float(r["z"]) for r in ds if r["set"] == s and r["z"] not in ("", "NA") and not r.get("site")}
        alt = [(g, zs[g]) for g in zs if g in resp and (sig_only or True)]
        altered = {r["gene"] for r in ds if r["set"] == s and r["z"] not in ("", "NA") and (sig_only or (r["p"] not in ("", "NA") and float(r["p"]) < alpha))}
        alt = [(g, z) for g, z in alt if g in altered]
        sz = statistics.pstdev(z for _, z in alt) if len(alt) > 1 else 1.0; sr = statistics.pstdev(resp[g] for g, _ in alt) if len(alt) > 1 else 1.0
        rows = []
        for g in slides.get(story, []):
            if g not in ok:
                continue
            z = zs.get(g) if g in altered else None; r = resp.get(g)
            score = -(z / (sz or 1.0)) * (r / (sr or 1.0)) if z is not None and r is not None else None
            rows.append({"story": story, "set": s, "node": g, "z": z, "response": r, "arm": arm, "tissue": tissue, "score": score, "links": red_blue(g)})
        rows.sort(key=lambda d: (d["score"] is None, -(d["score"] or 0)))
        picked += rows[:per_story]
    return picked


def super_list(out: Path, net: Network, per_story: int = 99, stories: Optional[Sequence[str]] = None) -> List[Dict[str, object]]:
    """The menu: one entry per node (best score first within story order), with all the stories it belongs to."""
    merged: Dict[str, Dict[str, object]] = {}
    for d in story_markers(out, net, per_story, stories):
        e = merged.setdefault(d["node"], {"node": d["node"], "stories": [], "best": d, "links": d["links"]})
        e["stories"].append(d)
        if (d["score"] or -1e9) > (e["best"]["score"] or -1e9):
            e["best"] = d
    return list(merged.values())


def describe(e: Dict[str, object]) -> str:
    """One menu line: the node, its slide(s), and — only when the data say so — the disease direction and whether the
    story's winning arm pushes it back."""
    parts = []
    for d in e["stories"]:
        tis = "muscle" if d["tissue"] == "muscle" else "blood"
        if d["z"] is None:
            parts.append(f"{d['story']} slide")
        else:
            way = "up" if d["z"] > 0 else "down"
            back = f", {ARM_NAME[d['arm']]} pushes it back in {tis}" if (d["score"] or 0) > 0 else ""
            parts.append(f"{d['story']}: {way} in disease (z {d['z']:+.0f}){back}")
    return f"{e['node']:<9} {e['links']} red/blue link{'s' if e['links'] != 1 else ''} · " + "; ".join(parts)


def write_menu(entries: List[Dict[str, object]], path: Path) -> None:
    """The menu as a CSV (one row per node x story), for the record."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["menu_no", "node", "story", "set", "disease_z", "arm", "tissue", "response", "reversal_score", "links"])
        for i, e in enumerate(entries, 1):
            for d in e["stories"]:
                rnd = lambda x, k: None if x is None else round(x, k)
                w.writerow([i, d["node"], d["story"], d["set"], rnd(d["z"], 3), d["arm"], d["tissue"], rnd(d["response"], 4), rnd(d["score"], 3), d["links"]])
