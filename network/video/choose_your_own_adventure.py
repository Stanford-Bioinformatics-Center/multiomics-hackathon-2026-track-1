#!/usr/bin/env python3
# =====================================================================================================
# video/choose_your_own_adventure.py — ONE PICTURE OF THE START-NODE CHOICES, TO POINT AT
# =====================================================================================================
#
# PURPOSE
#   The start-node menu of the music video as one simple picture titled "Choose Your Own Adventure": the nodes of the
#   blood future-T2D story slide (clean_for_slides/blood_t2d.png; the only slide with 4-node walks on figure 17 edges)
#   from which every 3-step walk stays on the slide and goes the full distance (exvideo/starts.py). Each choice is a
#   coloured dot with its first-step neighbours ON THE SLIDE in grey around it, so anyone can point at a node and
#   send its name. CYOA_SCOPE=<story> draws another slide's menu (as make_music_video.py --scope).
#
# HOW TO RUN
#   python3 network/video/choose_your_own_adventure.py            # -> $HACK_OUT/video/choose_your_own_adventure.png
#   Needs: the pipeline outputs (steps 01, 14, 20), R with ggplot2.
#
# OUTPUTS ($HACK_OUT/video/): start_menu.csv, start_menu_edges.csv, choose_your_own_adventure.png
# =====================================================================================================
from __future__ import annotations

import csv
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from exvideo import network, starts  # noqa: E402

GROUP = {("T2D muscle",): "T2D muscle", ("T2D blood",): "T2D blood", ("ageing blood",): "Ageing blood",
         ("T2D blood", "ageing blood"): "T2D + ageing blood"}


def main() -> int:
    out = Path(os.environ.get("HACK_OUT", str(Path.home() / "Desktop/output/hackathon-2026-track1/network")))
    net = network.load_network(out)
    scope = os.environ.get("CYOA_SCOPE", "T2D blood")                # the slide the walks stay on (as make_music_video --scope)
    menu = starts.super_list(out, net, stories=[scope])
    wnet = network.subnetwork(net, starts.scope_nodes(scope, net))    # first steps drawn = the ones the walk can take
    dest = out / "video"; dest.mkdir(parents=True, exist_ok=True)
    starts.write_menu(menu, dest / "start_menu.csv")
    with (dest / "start_menu_edges.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["choice", "group", "neighbour"])
        for e in menu:
            g = GROUP[tuple(sorted({d["story"] for d in e["stories"]}, key=["T2D muscle", "T2D blood", "ageing blood"].index))]
            for nb in sorted(x.other(e["node"]) for x in wnet.neighbours(e["node"])):
                w.writerow([e["node"], g, nb])
    png = dest / "choose_your_own_adventure.png"
    res = subprocess.run(["Rscript", str(HERE / "choose_your_own_adventure.R"), str(dest / "start_menu_edges.csv"), str(png)], capture_output=True, text=True)
    if res.returncode != 0:
        print(res.stderr[-800:], file=sys.stderr); return 1
    print(f"{len(menu)} choices -> {png}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
