"""The start-node menu: the most interesting markers of the analysis's three stories, as walk starts.

The three stories are step 20's primary results (20_stories.csv):
  T2D muscle    Öhman 2021 T2D muscle proteins; ENDURANCE reverses them (the endurance arm's muscle response)
  T2D blood     UK Biobank incident-T2D plasma proteins (Gadd 2024); RESISTANCE reverses them (blood protein response)
  ageing blood  UK Biobank ageing plasma proteins (Sun 2023); RESISTANCE reverses them (blood protein response)
For each story, exactly as step 20 does: the set's "altered" proteins (its own p threshold, 20_disease_sets.csv), the
disease z (20_disease_scores.csv), and the exercise response = the mean of the story tissue's cells in the winning
arm's node vectors (01_nodes_EE / RE.csv). A marker scores high when its disease change is large AND the winning arm
moves it the OTHER way: score = -(z / sd z) x (response / sd response) over the set's altered proteins (> 0 =
reversed). Only ROBUST starts are offered: figure 17 nodes from which every 3-step walk (never revisiting a node)
goes the full distance, whichever way it turns (walk.robust_starts). The super list merges the stories: a
node in several stories is listed once, with every story it belongs to.
"""
from __future__ import annotations

import csv
import statistics
from pathlib import Path
from typing import Dict, List

from .errors import InputError
from .network import Network
from .walk import robust_starts

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


def story_markers(out: Path, net: Network, per_story: int = 6) -> List[Dict[str, object]]:
    """The top `per_story` walkable markers of each story (ranked by reversal score)."""
    ok = set(robust_starts(net))
    meta = {r["set"]: r for r in _rows(out / "20_disease_sets.csv")}
    ds = _rows(out / "20_disease_scores.csv")
    resp_rows = {a: _rows(out / f"01_nodes_{a}.csv") for a in ("EE", "RE")}
    # red/blue links: the arm-specific edges the video draws bold (the figure 17 rule, as in network.build_facts)
    pool = sorted(abs(w) for e in net.edges.values() for w in (e.w_ee, e.w_re)); tau = pool[int(0.75 * (len(pool) - 1))] if pool else 0.0
    red_blue = lambda n: sum(1 for e in net.neighbours(n) if (e.w_ee >= tau) != (e.w_re >= tau))
    picked = []
    for story, s, tissue, arm in STORIES:
        if s not in meta:
            raise InputError(f"{s} is not in 20_disease_sets.csv (re-run step 20)")
        alpha, sig_only = float(meta[s]["alpha"]), meta[s]["sig_only"] == "TRUE"
        resp = _tissue_mean(resp_rows[arm], tissue)
        alt = [(r["gene"], float(r["z"])) for r in ds if r["set"] == s and r["z"] not in ("", "NA")
               and (sig_only or (r["p"] not in ("", "NA") and float(r["p"]) < alpha))]
        alt = [(g, z) for g, z in alt if g in resp]
        if len(alt) < 3:
            continue
        sz = statistics.pstdev(z for _, z in alt) or 1.0; sr = statistics.pstdev(resp[g] for g, _ in alt) or 1.0
        ranked = sorted(({"story": story, "set": s, "node": g, "z": z, "response": resp[g], "arm": arm, "tissue": tissue,
                          "score": -(z / sz) * (resp[g] / sr), "links": red_blue(g)}
                         for g, z in alt if g in ok), key=lambda d: -d["score"])
        picked += [d for d in ranked if d["score"] > 0][:per_story]
    return picked


def super_list(out: Path, net: Network, per_story: int = 6) -> List[Dict[str, object]]:
    """The menu: one entry per node (best score first within story order), with all the stories it belongs to."""
    merged: Dict[str, Dict[str, object]] = {}
    for d in story_markers(out, net, per_story):
        e = merged.setdefault(d["node"], {"node": d["node"], "stories": [], "best": d, "links": d["links"]})
        e["stories"].append(d)
        if d["score"] > e["best"]["score"]:
            e["best"] = d
    return list(merged.values())


def describe(e: Dict[str, object]) -> str:
    """One menu line: the node, its stories, and why (disease direction, the reversing arm's response)."""
    parts = []
    for d in e["stories"]:
        way = "up" if d["z"] > 0 else "down"
        tis = "muscle" if d["tissue"] == "muscle" else "blood"
        parts.append(f"{d['story']}: {way} in disease (z {d['z']:+.0f}), {ARM_NAME[d['arm']]} pushes it back in {tis}")
    return f"{e['node']:<9} {e['links']} red/blue link{'s' if e['links'] != 1 else ''} · " + "; ".join(parts)


def write_menu(entries: List[Dict[str, object]], path: Path) -> None:
    """The menu as a CSV (one row per node x story), for the record."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["menu_no", "node", "story", "set", "disease_z", "arm", "tissue", "response", "reversal_score", "links"])
        for i, e in enumerate(entries, 1):
            for d in e["stories"]:
                w.writerow([i, d["node"], d["story"], d["set"], round(d["z"], 3), d["arm"], d["tissue"], round(d["response"], 4), round(d["score"], 3), d["links"]])
