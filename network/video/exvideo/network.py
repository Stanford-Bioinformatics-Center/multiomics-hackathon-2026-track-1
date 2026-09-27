"""The joint network (step 14) and the figure 17 graph facts of the nodes and edges on a walk."""
from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, FrozenSet, List, Optional, Sequence

from .errors import InputError

ALPHA = 0.05            # adj. p threshold for "a significant exercise response", as in the figure 17 pages
ARM_SPECIFIC_Q = 0.75   # arm-specific edges: strong = |w| in the top 25% over both arms (figure 17 default)
DISEASE_SETS = {"ohman_2021": "T2D muscle proteome (Öhman 2021; z < 0 = lower in T2D)",
                "gadd_2024_ukb_incident_T2D": "UK Biobank plasma, future T2D (Gadd 2024; z > 0 = higher risk)",
                "sun_2023_ukb_age": "UK Biobank plasma, age (Sun 2023; z > 0 = higher with age)",
                "kjaergaard_2025_prot_pooled": "T2D muscle, Kjærgaard 2025 pooled (post hoc)"}


@dataclass(frozen=True)
class Edge:
    """One joint-network edge with its two per-arm weights."""
    a: str
    b: str
    edge_type: str
    w_ee: float
    w_re: float

    def other(self, node: str) -> str:
        """The node at the other end of this edge."""
        return self.b if node == self.a else self.a


@dataclass
class Network:
    """The joint network: node table rows and edges keyed by their (unordered) node pair."""
    nodes: Dict[str, Dict[str, str]]
    edges: Dict[FrozenSet[str], Edge]

    def neighbours(self, node: str) -> List[Edge]:
        """Every edge touching `node`, in a fixed order (sorted by the other end's name)."""
        return sorted((e for k, e in self.edges.items() if node in k), key=lambda e: e.other(node))


def read_csv(path: Path) -> List[Dict[str, str]]:
    """Read a CSV into row dictionaries; stop with an InputError if it is missing."""
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
    """Round to 3 significant digits for the prompt (None stays None)."""
    return None if x is None else float(f"{x:.3g}")


def load_network(out: Path) -> Network:
    """The joint network from step 14 (14_joint_nodes.csv, 14_joint_edges.csv)."""
    edges = {frozenset((r["node_a"], r["node_b"])): Edge(r["node_a"], r["node_b"], r["edge_type"], float(r["w_EE"]), float(r["w_RE"]))
             for r in read_csv(out / "14_joint_edges.csv")}
    nodes = {r["node"]: r for r in read_csv(out / "14_joint_nodes.csv")}
    return Network(nodes, edges)


def subnetwork(net: Network, nodes) -> Network:
    """The network restricted to `nodes` (the edges between them only), e.g. one story slide."""
    keep = {n for n in nodes if n in net.nodes}
    return Network({n: net.nodes[n] for n in keep}, {k: e for k, e in net.edges.items() if k <= keep})


def build_facts(walk: Sequence[str], net: Network, out: Path) -> Dict[str, Any]:
    """What the figure 17 data say about every node and edge of the walk (the walk must already be checked)."""
    # Arm-specific threshold over ALL joint edges and both arms (the figure 17 rule), and the hub ranking.
    pool = sorted(abs(w) for e in net.edges.values() for w in (e.w_ee, e.w_re))
    tau = pool[int(ARM_SPECIFIC_Q * (len(pool) - 1))]
    ranked = sorted(net.nodes, key=lambda n: -max(float(net.nodes[n]["strength_EE"]), float(net.nodes[n]["strength_RE"])))
    string_score = {frozenset((r["symbol_a"], r["symbol_b"])): num(r["combined_score"]) for r in read_csv(out / "02_edges.csv")}
    module = {r["node"]: r["module"] for r in read_csv(out / "17_modules.csv") if r["network"] == "joint"}
    mname = {r["module"]: r for r in read_csv(out / "17_module_names.csv") if r["network"] == "joint"}
    cells = [r for r in read_csv(out / "17_node_cell_stats.csv") if r["node"] in walk and r["arm"] in ("EE", "RE")]
    phos = [r for r in read_csv(out / "17_phospho_site_stats.csv") if r["protein"] in walk and r["arm"] in ("EE", "RE")]
    gly = {r["protein"]: r for r in read_csv(out / "17_glygen_protein_annotation.csv") if r["protein"] in walk}
    dz = {(r["set"], r["gene"]): r for r in read_csv(out / "20_disease_scores.csv") if r["set"] in DISEASE_SETS and r["gene"] in walk and not r["site"]}
    # What the VIDEO SHOWS of a node's links (the figure 17 page in arm-specific mode, render_walk.js): only the
    # arm-specific edges are drawn (bold red = strong after endurance only, blue = resistance only; every other edge is a
    # faint thin grey line), plus the walk's own edges in gold. Link COUNTS are the red + blue links only (Vidal,
    # 2026-09-27); the gold walk links are listed for the hand-offs but not counted. The full joint-network degree is not
    # given to the model, so the lyrics never count links the viewer cannot see.
    walk_pairs = {frozenset(p) for p in zip(walk, walk[1:])}

    def arm_of(e: Edge) -> str:
        hE, hR = e.w_ee >= tau, e.w_re >= tau
        return "both" if hE and hR else "EE" if hE else "RE" if hR else "neither"

    def on_screen(n: str) -> Dict[str, Any]:
        es = net.neighbours(n)
        ee = sorted((e for e in es if arm_of(e) == "EE"), key=lambda e: -e.w_ee)
        re_ = sorted((e for e in es if arm_of(e) == "RE"), key=lambda e: -e.w_re)
        walked = [e.other(n) for e in es if frozenset((n, e.other(n))) in walk_pairs]
        return {"red_blue_count": len(ee) + len(re_), "endurance_only_links_red": [e.other(n) for e in ee][:8],
                "resistance_only_links_blue": [e.other(n) for e in re_][:8], "walk_links_gold_not_counted": walked}

    def best_cell(node: str, arm: str) -> Optional[Dict[str, Any]]:
        """The node's most significant exercise response in one arm (adj. p < ALPHA), or None."""
        sig = [c for c in cells if c["node"] == node and c["arm"] == arm and (num(c["adj_p"]) or 1.0) < ALPHA]
        if not sig:
            return None
        c = min(sig, key=lambda c: num(c["adj_p"]) or 1.0)
        return {"tissue": c["tissue"], "ome": c["ome"], "time": c["time"], "logFC": r3(num(c["logFC"])), "adj_p": r3(num(c["adj_p"])), "n_significant_cells": len(sig)}

    def phospho(node: str) -> Dict[str, Any]:
        """Measured and responding MoTrPAC muscle phosphosites of a protein, by arm."""
        by_site: Dict[str, set] = {}
        for r in phos:
            if r["protein"] == node and r["tissue"] == "muscle" and (num(r["adj_p"]) or 1.0) < ALPHA:
                by_site.setdefault(r["site"], set()).add(r["arm"])
        cls = {s: ("both" if len(a) == 2 else ("endurance only" if "EE" in a else "resistance only")) for s, a in by_site.items()}
        return {"measured_sites": len({r["site"] for r in phos if r["protein"] == node}),
                "responding_sites": {k: sorted(s for s, c in cls.items() if c == k)[:6] for k in ("endurance only", "resistance only", "both")}}

    node_facts: List[Dict[str, Any]] = []
    for n in walk:
        nd = net.nodes[n]; m = module.get(n); mn = mname.get(m or "", {}); g = gly.get(n, {})
        node_facts.append({
            "node": n, "type": nd["node_type"], "class": nd.get("class") or None,
            "links_on_screen": on_screen(n),
            "strength_endurance": r3(float(nd["strength_EE"])), "strength_resistance": r3(float(nd["strength_RE"])),
            "hub_rank_in_network": ranked.index(n) + 1, "network_size": len(net.nodes),
            "mean_response_endurance": r3(num(nd["resp_EE"])), "mean_response_resistance": r3(num(nd["resp_RE"])),
            "module": m, "module_pathway": mn.get("path_name") or None, "module_pathway_fdr": r3(num(mn.get("path_fdr"))),
            "strongest_significant_response_endurance": best_cell(n, "EE"),
            "strongest_significant_response_resistance": best_cell(n, "RE"),
            "disease_and_ageing": {label: {"z": r3(num(dz[(s, n)]["z"])), "p": r3(num(dz[(s, n)]["p"]))} for s, label in DISEASE_SETS.items() if (s, n) in dz},
            "phosphosites_muscle": phospho(n) if nd["node_type"] == "protein" else None,
            "glycosylation_sites": {"N_linked": g.get("glyco_N_sites"), "O_linked": g.get("glyco_O_sites"), "O_GlcNAc": g.get("glyco_OGlcNAc_sites")} if g else None,
            "kinase": {"is_kinase": g.get("is_kinase") == "1", "substrate_sites": g.get("substrate_sites"), "known_phosphosites": g.get("phosphosites")} if g else None,
        })
    edge_facts: List[Dict[str, Any]] = []
    for u, v in zip(walk, walk[1:]):
        e = net.edges[frozenset((u, v))]
        spec = ("both arms" if e.w_ee >= tau and e.w_re >= tau else "endurance only" if e.w_ee >= tau else "resistance only" if e.w_re >= tau else "neither arm")
        edge_facts.append({"from": u, "to": v, "type": e.edge_type, "physical_evidence": {"STRING_combined_score": string_score.get(frozenset((u, v)))},
                           "weight_endurance": r3(e.w_ee), "weight_resistance": r3(e.w_re), "weight_difference_EE_minus_RE": r3(e.w_ee - e.w_re),
                           "strong_in": spec, "strong_threshold_abs_w": r3(tau)})
    return {"walk": list(walk),
            "links_note": ("links_on_screen is what the video shows. red_blue_count counts the RED (strong after endurance only) and "
                           "BLUE (strong after resistance only) links, the only links drawn bold; when a line counts a node's links, "
                           "use red_blue_count and nothing else. The walk's own links (gold) are for the hand-off to the next node, "
                           "not for counting. The network has other, faint links that are not shown: never state a node's total "
                           "number of links or partners."),
            "about": ("Joint network of the MoTrPAC acute-exercise response (471 genes/proteins, 450 metabolites; adipose, blood, "
                      "muscle). An edge exists only if STRING (score >= 700) or Rhea physically links the two molecules; its weight "
                      "per arm is the dot product of their normalised exercise responses (positive = they respond together)."),
            "nodes": node_facts, "edges": edge_facts}


def fact_card(node_fact: Dict[str, Any]) -> List[str]:
    """Three short on-screen lines about a node, from its facts (for the video's fact card)."""
    f = node_fact; k = f["links_on_screen"]["red_blue_count"]
    lines = [f"hub rank {f['hub_rank_in_network']} of {f['network_size']} · {k} red/blue link{'s' if k != 1 else ''}"]
    lines.append(f"strength: endurance {f['strength_endurance']} · resistance {f['strength_resistance']}")
    d = f.get("disease_and_ageing") or {}
    ukb = [(k.split(",")[1].split("(")[0].strip(), v["z"]) for k, v in d.items() if k.startswith("UK Biobank")]
    if ukb:
        lines.append("UK Biobank " + " · ".join(f"{k} z {z:+.1f}" for k, z in ukb))
    elif f.get("module_pathway"):
        lines.append(f"module: {f['module_pathway']}")
    return lines
