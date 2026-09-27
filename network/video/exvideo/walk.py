"""Walks through the joint network: a fixed walk (checked) or a random walk from a prompted node."""
from __future__ import annotations

import csv
import difflib
import os
import random
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

from .errors import WalkError
from .network import Network


def resolve_node(name: str, net: Network) -> str:
    """Match a typed node name to the network (exact, then case-insensitive); suggest close names if none."""
    name = name.strip()
    if name in net.nodes:
        return name
    folded = {n.lower(): n for n in net.nodes}
    if name.lower() in folded:
        return folded[name.lower()]
    close = difflib.get_close_matches(name, list(net.nodes), n=5, cutoff=0.6) or difflib.get_close_matches(name.lower(), list(folded), n=5, cutoff=0.6)
    hint = f"; did you mean: {', '.join(folded.get(c, c) for c in close)}?" if close else ""
    raise WalkError(f"{name!r} is not a node of the joint network (353 nodes with at least one edge){hint}")


def check_walk(walk: Sequence[str], net: Network) -> None:
    """Every node must be in the network, none twice, and every consecutive pair must be a (physical) edge."""
    if len(walk) < 2:
        raise WalkError("a walk needs at least two nodes")
    if len(set(walk)) != len(walk):
        raise WalkError(f"the walk visits a node twice: {' -> '.join(walk)}")
    unknown = [n for n in walk if n not in net.nodes]
    if unknown:
        raise WalkError(f"not in the joint network: {', '.join(unknown)}")
    for u, v in zip(walk, walk[1:]):
        if frozenset((u, v)) not in net.edges:
            raise WalkError(f"{u} - {v} is not an edge of the joint network (no STRING / Rhea link), so the walk cannot take it")


TEAM_WALKER = Path(__file__).resolve().parents[3] / "random_walk" / "random_walks.R"


def team_walk(start: str, arm: str, seed: Optional[int], out: Path, work: Path) -> Tuple[Tuple[str, ...], int, List[float]]:
    """A 3-step walk from `start` with the TEAM's walker (random_walk/random_walks.R, Rscript).

    Its rule: each step goes to a not-yet-visited neighbour with probability set by the chosen arm's edge weights
    (EE = endurance, RE = resistance); dead ends are redrawn (see random_walk/README.md). The walk is then checked here
    against the physical edges. Returns (walk, seed, step probabilities).
    """
    if not TEAM_WALKER.exists():
        raise WalkError(f"the team's walker is not in the repository ({TEAM_WALKER}); use --walker builtin")
    if shutil.which("Rscript") is None:
        raise WalkError("Rscript not found (needed by the team's walker); use --walker builtin")
    seed = random.SystemRandom().randrange(1, 10**6) if seed is None else seed
    work.mkdir(parents=True, exist_ok=True)
    res = subprocess.run(["Rscript", str(TEAM_WALKER), start, arm, str(seed)], capture_output=True, text=True,
                         env={**os.environ, "HACK_OUT": str(out), "RW_OUT": str(work)})
    f = work / "18_walk.csv"
    if res.returncode != 0 or not f.exists():
        raise WalkError(f"the team's walker failed: {(res.stderr or res.stdout).strip()[-400:]}")
    with f.open(newline="", encoding="utf-8") as fh:
        rows = sorted(csv.DictReader(fh), key=lambda r: int(r["step"]))
    return tuple(r["node"] for r in rows), seed, [float(r["p_step"]) for r in rows]


def random_walk(start: str, steps: int, net: Network, seed: Optional[int] = None, tries: int = 200) -> Tuple[Tuple[str, ...], int]:
    """A random walk of `steps` edges from `start` along the network's physical edges.

    Rule: at each step, move to a not-yet-visited neighbour with probability proportional to the edge's strength,
    max(|w_EE|, |w_RE|) — strongly co-regulated links are likelier. A walk that reaches a dead end (no unvisited
    neighbour) is discarded and redrawn, up to `tries` times. Reproducible: the same `seed` gives the same walk.
    Returns (walk, seed used).
    """
    if steps < 1:
        raise WalkError("a random walk needs at least one step")
    if not net.neighbours(start):
        raise WalkError(f"{start} has no edges, so no walk can start there")
    seed = random.SystemRandom().randrange(1, 10**6) if seed is None else seed
    rng = random.Random(seed)
    for _ in range(tries):
        walk: List[str] = [start]
        while len(walk) <= steps:
            options = [e for e in net.neighbours(walk[-1]) if e.other(walk[-1]) not in walk]
            if not options:
                break
            weights = [max(abs(e.w_ee), abs(e.w_re)) + 1e-9 for e in options]
            walk.append(rng.choices(options, weights=weights, k=1)[0].other(walk[-1]))
        if len(walk) == steps + 1:
            check_walk(walk, net)
            return tuple(walk), seed
    raise WalkError(f"no {steps}-step walk without revisits found from {start} in {tries} tries (its neighbourhood is too small)")
