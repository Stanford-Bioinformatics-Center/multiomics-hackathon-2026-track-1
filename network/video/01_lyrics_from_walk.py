#!/usr/bin/env python3
# =====================================================================================================
# video/01_lyrics_from_walk.py — MUSIC VIDEO STAGE 1: FROM A WALK THROUGH THE NETWORK TO 16 BARS OF LYRICS
# =====================================================================================================
#
# PURPOSE
#   The first, "upstream" stage of the Team 2-PAC music video. It takes a walk through our joint network (for now a
#   fixed four-node walk, HYOU1 -> HSP90B1 -> CDC37 -> SRC), gathers what the figure 17 graph data say about every
#   node and every edge on it, and asks Claude — with the team's own prompt, word for word — for 16 bars of lyrics,
#   4 per node, in walk order. The lyrics are saved in a structured file that the video stage (the teammate's walk /
#   video script) reads to know which node to show while each bar plays.
#
# WHAT THIS SCRIPT DOES (plain language)
#   1. The walk: taken from --walk (default HYOU1,HSP90B1,CDC37,SRC). Every consecutive pair must be an edge of the
#      joint network (a physical STRING / Rhea link — the network's hard gate); otherwise the script stops, because
#      the video would show a connection that does not exist. (--start / --steps are reserved for the planned random
#      walk from a prompted node; not implemented yet, see random_walk().)
#   2. The facts (the same tables the figure 17 pages are built from, step 14 / 17 / 20 outputs):
#        per node:  type, degree, strength per arm (sum of |edge weight|) and its hub rank, network module and the
#                   module's pathway name, mean response per arm, the strongest significant exercise response per arm
#                   (tissue, ome, time, logFC, adj. p), T2D and ageing directions (Öhman 2021 muscle T2D, UK Biobank
#                   future T2D and age, Kjærgaard 2025 pooled), responding MoTrPAC phosphosites per arm, glycosylation;
#        per edge:  type, STRING score, weight after endurance and after resistance, their difference, and whether
#                   it is an arm-specific edge (strong = top 25% of |w| over all joint edges and both arms: the rule of
#                   the figure 17 "arm-specific edges" switch and exnet::arm_specific_edges()).
#   3. The prompt: the team's prompt (PROMPT_TEMPLATE below, verbatim, with the walk filled in) + the facts as JSON +
#      the output format (JSON: 16 bars, each tagged with its node).
#   4. Claude: --backend cli (default) runs the Claude Code command line (`claude -p`), --backend api calls the
#      Anthropic Messages API (needs ANTHROPIC_API_KEY), --backend none only writes the prompt. The answer is checked
#      (16 bars, 4 per node, in walk order); if it is malformed, the script asks once more, saying what was wrong.
#
# HOW TO RUN
#   python3 network/video/01_lyrics_from_walk.py                         # default walk, Claude Code CLI
#   python3 network/video/01_lyrics_from_walk.py --backend none          # write the prompt only (no model call)
#   python3 network/video/01_lyrics_from_walk.py --walk HYOU1,HSP90B1,CDC37,SRC --model claude-opus-5-5
#   Needs the pipeline outputs (bash network/run_all.sh) and, for --backend cli, Claude Code on the PATH.
#
# INPUTS:  $HACK_OUT/14_joint_edges.csv, 14_joint_nodes.csv, 02_edges.csv, 17_modules.csv, 17_module_names.csv,
#          17_node_cell_stats.csv, 17_phospho_site_stats.csv, 17_glygen_protein_annotation.csv, 20_disease_scores.csv
# OUTPUTS: $HACK_OUT/video/<walk>/walk_facts.json   the facts sent to the model
#          $HACK_OUT/video/<walk>/lyrics_prompt.md  the exact prompt
#          $HACK_OUT/video/<walk>/lyrics_raw.txt    the model's raw answer
#          $HACK_OUT/video/<walk>/lyrics.json       CONTRACT for the video stage: {"walk", "title", "bars": [{"bar",
#                                                   "node", "text"}], "model", "backend", "created"}
#          $HACK_OUT/video/<walk>/lyrics.md         the lyrics, readable
#
# KNOWN LIMITS
#   A language model writes the lyrics, so two runs give different lyrics (the facts and the prompt are fixed and
#   saved, so every run is traceable). The lyrics are a creative summary of the data, not a scientific claim.
# =====================================================================================================
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import urllib.error
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

# ---- the team's prompt (verbatim; {walk} is replaced by the walk, e.g. "HYOU1 -> HSP90B1 -> CDC37 -> SRC") ----------
PROMPT_TEMPLATE = (
    "give me 16 bars of 2pac rap lyrics which summarize the most interesting story from these 4 nodes: {walk}. "
    "4 bars per node, and write the lyrics based on what is most interesting based on the figure 17 graph data. "
    "right now, we will keep the walk listed ({walk}) and later we will extend the script to do a random walk on a "
    "prompted node"
)
# What is appended after the team's prompt: the data, originality, and a format the video stage can read.
PROMPT_SUFFIX = """

Write ORIGINAL lyrics in the style of 1990s West Coast rap (do not quote or reuse lines from any existing song).
Use the figure 17 graph data below: pick what is most interesting (strong or arm-specific edges, the hub, disease and
ageing directions, phosphosites that respond to one kind of exercise) and keep every biological statement true to the
data. Bars 1-4 are about node 1, bars 5-8 node 2, bars 9-12 node 3, bars 13-16 node 4; mention the link to the next
node where it fits.

FIGURE 17 GRAPH DATA (JSON):
{facts}

Answer with JSON only, exactly this shape (no other text):
{{"title": "<song title>", "bars": [{{"bar": 1, "node": "<node>", "text": "<one bar>"}}, ... 16 items ...]}}
"""
DEFAULT_WALK = ("HYOU1", "HSP90B1", "CDC37", "SRC")
BARS_PER_NODE = 4
DEFAULT_MODEL = "claude-opus-5-5"
ALPHA = 0.05            # adj. p threshold for "a significant exercise response", as in the figure 17 pages
ARM_SPECIFIC_Q = 0.75   # arm-specific edges: strong = |w| in the top 25% over both arms (figure 17 default)


# ---- errors: every failure says what is wrong ------------------------------------------------------------------------
class VideoStageError(Exception):
    """Base class for problems this stage detects."""


class WalkError(VideoStageError):
    """The walk is not a valid path through the joint network."""


class InputError(VideoStageError):
    """A pipeline output this stage needs is missing or malformed."""


class ModelError(VideoStageError):
    """The model could not be called, or did not return 16 valid bars."""


# ---- small typed helpers ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Edge:
    """One joint-network edge with its two per-arm weights."""
    a: str
    b: str
    edge_type: str
    w_ee: float
    w_re: float


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Read a CSV into a list of row dictionaries; stop with an InputError if it is missing."""
    if not path.exists():
        raise InputError(f"missing pipeline output {path} (run: bash network/run_all.sh)")
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def num(x: Optional[str]) -> Optional[float]:
    """Text to a finite number, or None (empty, NA, not a number)."""
    try:
        v = float(x) if x not in (None, "", "NA") else None
    except ValueError:
        return None
    return v if v is not None and math.isfinite(v) else None


def r3(x: Optional[float]) -> Optional[float]:
    """Round to 3 significant decimals for the prompt (None stays None)."""
    return None if x is None else float(f"{x:.3g}")


# ---- the walk --------------------------------------------------------------------------------------------------------
def random_walk(start: str, steps: int, seed: int) -> Tuple[str, ...]:
    """PLANNED (not implemented yet): a random walk of `steps` edges from the prompted node `start`.

    Intended rule: at each step move to a neighbour with probability proportional to |w| of the connecting edge
    (optionally in one arm), never revisiting a node, reproducible with `seed`; the result is passed to the same
    checks and facts as a fixed walk.
    """
    raise NotImplementedError("random walks from a prompted node are planned; use --walk for now")


def check_walk(walk: Sequence[str], edges: Dict[frozenset, Edge], nodes: Dict[str, Dict[str, str]]) -> List[Edge]:
    """Every node must be in the joint network and every consecutive pair must be one of its (physical) edges."""
    if len(walk) < 2:
        raise WalkError("a walk needs at least two nodes")
    if len(set(walk)) != len(walk):
        raise WalkError(f"the walk visits a node twice: {' -> '.join(walk)}")
    unknown = [n for n in walk if n not in nodes]
    if unknown:
        raise WalkError(f"not in the joint network: {', '.join(unknown)}")
    steps: List[Edge] = []
    for u, v in zip(walk, walk[1:]):
        e = edges.get(frozenset((u, v)))
        if e is None:
            raise WalkError(f"{u} - {v} is not an edge of the joint network (no STRING / Rhea link), so the walk cannot take it")
        steps.append(e)
    return steps


# ---- the facts (figure 17 graph data) --------------------------------------------------------------------------------
def build_facts(walk: Sequence[str], out: Path) -> Dict[str, Any]:
    """Collect, for every node and edge of the walk, what the figure 17 data say (see the header)."""
    joint = read_csv(out / "14_joint_edges.csv")
    edges = {frozenset((r["node_a"], r["node_b"])): Edge(r["node_a"], r["node_b"], r["edge_type"], float(r["w_EE"]), float(r["w_RE"])) for r in joint}
    nodes = {r["node"]: r for r in read_csv(out / "14_joint_nodes.csv")}
    steps = check_walk(walk, edges, nodes)

    # Arm-specific threshold over ALL joint edges and both arms (the figure 17 rule).
    pool = sorted(abs(w) for e in edges.values() for w in (e.w_ee, e.w_re))
    tau = pool[int(ARM_SPECIFIC_Q * (len(pool) - 1))]
    # Hub rank: nodes ordered by their larger per-arm strength.
    ranked = sorted(nodes, key=lambda n: -max(float(nodes[n]["strength_EE"]), float(nodes[n]["strength_RE"])))
    # STRING scores by gene symbol pair.
    string_score = {frozenset((r["symbol_a"], r["symbol_b"])): num(r["combined_score"]) for r in read_csv(out / "02_edges.csv")}
    # Modules (joint network) and their names.
    module = {r["node"]: r["module"] for r in read_csv(out / "17_modules.csv") if r["network"] == "joint"}
    mname = {r["module"]: r for r in read_csv(out / "17_module_names.csv") if r["network"] == "joint"}
    # Per-cell statistics, phosphosites, glycosylation, disease / ageing directions.
    cells = [r for r in read_csv(out / "17_node_cell_stats.csv") if r["node"] in walk and r["arm"] in ("EE", "RE")]
    phos = [r for r in read_csv(out / "17_phospho_site_stats.csv") if r["protein"] in walk and r["arm"] in ("EE", "RE")]
    gly = {r["protein"]: r for r in read_csv(out / "17_glygen_protein_annotation.csv") if r["protein"] in walk}
    sets = {"ohman_2021": "T2D muscle proteome (Öhman 2021; z < 0 = lower in T2D)",
            "gadd_2024_ukb_incident_T2D": "UK Biobank plasma, future T2D (Gadd 2024; z > 0 = higher risk)",
            "sun_2023_ukb_age": "UK Biobank plasma, age (Sun 2023; z > 0 = higher with age)",
            "kjaergaard_2025_prot_pooled": "T2D muscle, Kjærgaard 2025 pooled (post hoc)"}
    dz = {(r["set"], r["gene"]): r for r in read_csv(out / "20_disease_scores.csv") if r["set"] in sets and r["gene"] in walk and not r["site"]}

    def best_cell(node: str, arm: str) -> Optional[Dict[str, Any]]:
        """The node's most significant exercise response in one arm (adj. p < ALPHA), or None."""
        sig = [c for c in cells if c["node"] == node and c["arm"] == arm and (num(c["adj_p"]) or 1.0) < ALPHA]
        if not sig:
            return None
        c = min(sig, key=lambda c: num(c["adj_p"]) or 1.0)
        return {"tissue": c["tissue"], "ome": c["ome"], "time": c["time"], "logFC": r3(num(c["logFC"])), "adj_p": r3(num(c["adj_p"])),
                "n_significant_cells": len(sig)}

    def phospho(node: str) -> Dict[str, Any]:
        """Responding MoTrPAC muscle phosphosites of a protein: which arm, and the strongest few."""
        by_site: Dict[str, set] = {}
        for r in phos:
            if r["protein"] == node and r["tissue"] == "muscle" and (num(r["adj_p"]) or 1.0) < ALPHA:
                by_site.setdefault(r["site"], set()).add(r["arm"])
        cls = {s: ("both" if len(a) == 2 else ("endurance only" if "EE" in a else "resistance only")) for s, a in by_site.items()}
        return {"measured_sites": len({r["site"] for r in phos if r["protein"] == node}),
                "responding_sites": {k: sorted(s for s, c in cls.items() if c == k)[:6] for k in ("endurance only", "resistance only", "both")}}

    node_facts = []
    for n in walk:
        nd = nodes[n]; m = module.get(n); mn = mname.get(m or "", {}); g = gly.get(n, {})
        node_facts.append({
            "node": n, "type": nd["node_type"], "degree": int(nd["degree"]),
            "strength_endurance": r3(float(nd["strength_EE"])), "strength_resistance": r3(float(nd["strength_RE"])),
            "hub_rank_in_network": ranked.index(n) + 1, "network_size": len(nodes),
            "mean_response_endurance": r3(num(nd["resp_EE"])), "mean_response_resistance": r3(num(nd["resp_RE"])),
            "module": m, "module_pathway": mn.get("path_name") or None, "module_pathway_fdr": r3(num(mn.get("path_fdr"))),
            "strongest_significant_response_endurance": best_cell(n, "EE"),
            "strongest_significant_response_resistance": best_cell(n, "RE"),
            "disease_and_ageing": {label: {"z": r3(num(dz[(s, n)]["z"])), "p": r3(num(dz[(s, n)]["p"]))} for s, label in sets.items() if (s, n) in dz},
            "phosphosites_muscle": phospho(n),
            "glycosylation_sites": {"N_linked": g.get("glyco_N_sites"), "O_linked": g.get("glyco_O_sites"), "O_GlcNAc": g.get("glyco_OGlcNAc_sites")} if g else None,
            "is_kinase": (g.get("is_kinase") == "1") if g else None,
        })
    edge_facts = []
    for e, (u, v) in zip(steps, zip(walk, walk[1:])):
        spec = ("both arms" if e.w_ee >= tau and e.w_re >= tau else "endurance only" if e.w_ee >= tau else "resistance only" if e.w_re >= tau else "neither arm")
        edge_facts.append({"from": u, "to": v, "type": e.edge_type, "physical_evidence": {"STRING_combined_score": string_score.get(frozenset((u, v)))},
                           "weight_endurance": r3(e.w_ee), "weight_resistance": r3(e.w_re), "weight_difference_EE_minus_RE": r3(e.w_ee - e.w_re),
                           "strong_in": spec, "strong_threshold_abs_w": r3(tau)})
    return {"walk": list(walk),
            "about": ("Joint network of the MoTrPAC acute-exercise response (471 genes/proteins, 450 metabolites; adipose, blood, "
                      "muscle). An edge exists only if STRING (score >= 700) or Rhea physically links the two molecules; its weight "
                      "per arm is the dot product of their normalised exercise responses (positive = they respond together)."),
            "nodes": node_facts, "edges": edge_facts}


# ---- the model call --------------------------------------------------------------------------------------------------
def call_cli(prompt: str, model: str) -> str:
    """Ask Claude through the Claude Code command line (non-interactive: claude -p)."""
    exe = shutil.which("claude")
    if exe is None:
        raise ModelError("Claude Code ('claude') is not on the PATH; install it or use --backend api / none")
    try:
        res = subprocess.run([exe, "-p", prompt, "--model", model, "--output-format", "text"], capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as err:
        raise ModelError("the Claude Code call timed out (600 s)") from err
    if res.returncode != 0:
        raise ModelError(f"claude -p failed (exit {res.returncode}): {res.stderr.strip()[:500]}")
    return res.stdout


def call_api(prompt: str, model: str) -> str:
    """Ask Claude through the Anthropic Messages API (needs ANTHROPIC_API_KEY)."""
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        raise ModelError("ANTHROPIC_API_KEY is not set; use --backend cli or none")
    body = json.dumps({"model": model, "max_tokens": 2000, "messages": [{"role": "user", "content": prompt}]}).encode()
    req = urllib.request.Request("https://api.anthropic.com/v1/messages", data=body, method="POST",
                                 headers={"x-api-key": key, "anthropic-version": "2023-06-01", "content-type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            out = json.loads(r.read().decode())
    except urllib.error.HTTPError as err:
        raise ModelError(f"Anthropic API error {err.code}: {err.read().decode()[:500]}") from err
    except urllib.error.URLError as err:
        raise ModelError(f"could not reach the Anthropic API: {err.reason}") from err
    return "".join(block.get("text", "") for block in out.get("content", []) if block.get("type") == "text")


def parse_bars(text: str, walk: Sequence[str]) -> Dict[str, Any]:
    """Pull the JSON out of the answer and check it: 16 bars, numbered 1-16, 4 per node in walk order, non-empty."""
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        raise ModelError("the answer contains no JSON object")
    try:
        obj = json.loads(m.group(0))
    except json.JSONDecodeError as err:
        raise ModelError(f"the answer's JSON does not parse: {err}") from err
    bars = obj.get("bars")
    need = BARS_PER_NODE * len(walk)
    if not isinstance(bars, list) or len(bars) != need:
        raise ModelError(f"expected {need} bars, got {len(bars) if isinstance(bars, list) else 'none'}")
    for i, b in enumerate(bars):
        want = walk[i // BARS_PER_NODE]
        if not isinstance(b, dict) or not str(b.get("text", "")).strip():
            raise ModelError(f"bar {i + 1} is empty or malformed")
        if b.get("node") != want:
            raise ModelError(f"bar {i + 1} is tagged {b.get('node')!r} but should be about {want!r} (4 bars per node, in walk order)")
        b["bar"] = i + 1
    return {"title": str(obj.get("title") or "Untitled"), "bars": [{"bar": b["bar"], "node": b["node"], "text": str(b["text"]).strip()} for b in bars]}


# ---- main --------------------------------------------------------------------------------------------------------------
def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build the facts and the prompt, ask the model, check and save the lyrics."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--walk", default=",".join(DEFAULT_WALK), help="comma-separated nodes (default: %(default)s)")
    ap.add_argument("--start", help="(planned) start node of a random walk")
    ap.add_argument("--steps", type=int, default=3, help="(planned) number of random-walk steps")
    ap.add_argument("--seed", type=int, default=20260926, help="(planned) random-walk seed")
    ap.add_argument("--backend", choices=("cli", "api", "none"), default="cli", help="how to ask Claude (default: %(default)s)")
    ap.add_argument("--model", default=DEFAULT_MODEL, help="Claude model (default: %(default)s)")
    ap.add_argument("--out", default=os.environ.get("HACK_OUT", str(Path.home() / "Desktop/output/hackathon-2026-track1/network")), help="pipeline output folder")
    args = ap.parse_args(argv)
    try:
        walk = random_walk(args.start, args.steps, args.seed) if args.start else tuple(s.strip() for s in args.walk.split(",") if s.strip())
        out = Path(args.out)
        facts = build_facts(walk, out)
        prompt = PROMPT_TEMPLATE.format(walk=" -> ".join(walk)) + PROMPT_SUFFIX.format(facts=json.dumps(facts, indent=2, ensure_ascii=False))
        dest = out / "video" / "_".join(walk); dest.mkdir(parents=True, exist_ok=True)
        (dest / "walk_facts.json").write_text(json.dumps(facts, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (dest / "lyrics_prompt.md").write_text(prompt + "\n", encoding="utf-8")
        print(f"walk {' -> '.join(walk)}: facts and prompt -> {dest}")
        if args.backend == "none":
            return 0
        ask = call_cli if args.backend == "cli" else call_api
        raw = ask(prompt, args.model)
        try:
            lyr = parse_bars(raw, walk)
        except ModelError as first:                      # one retry, telling the model what was wrong
            print(f"answer rejected ({first}); asking once more", file=sys.stderr)
            raw = ask(prompt + f"\n\nYour previous answer was rejected: {first}. Return only the JSON, exactly as specified.", args.model)
            lyr = parse_bars(raw, walk)
        (dest / "lyrics_raw.txt").write_text(raw, encoding="utf-8")
        rec = {"walk": list(walk), **lyr, "model": args.model, "backend": args.backend, "created": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        (dest / "lyrics.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        md = [f"# {lyr['title']}", "", f"*Walk: {' -> '.join(walk)} · {args.model}*", ""]
        for i, n in enumerate(walk):
            md += [f"**{n}**", ""] + [b["text"] for b in lyr["bars"][i * BARS_PER_NODE:(i + 1) * BARS_PER_NODE]] + [""]
        (dest / "lyrics.md").write_text("\n".join(md), encoding="utf-8")
        print(f"lyrics ({len(lyr['bars'])} bars) -> {dest / 'lyrics.json'}")
        return 0
    except VideoStageError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
