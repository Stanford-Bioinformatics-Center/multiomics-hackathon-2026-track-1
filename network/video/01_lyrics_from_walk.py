#!/usr/bin/env python3
# =====================================================================================================
# video/01_lyrics_from_walk.py — MUSIC VIDEO, LYRICS ONLY: walk -> figure 17 facts -> the team's prompt -> 16 bars
# =====================================================================================================
# The first stage of make_music_video.py on its own (no Suno, no video): useful to try walks and lyrics quickly.
# HOW TO RUN
#   python3 network/video/01_lyrics_from_walk.py                                  # the default walk
#   python3 network/video/01_lyrics_from_walk.py --start HYOU1 --seed 7           # a random walk from HYOU1
#   python3 network/video/01_lyrics_from_walk.py --backend none                   # facts + prompt only
# OUTPUTS: $HACK_OUT/video/<walk>/walk_facts.json, lyrics_prompt.md, lyrics_raw.txt, lyrics.json, lyrics.md,
#          suno_lyrics.txt, suno_style.txt (see video/README.md)
# =====================================================================================================
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exvideo import lyrics, network, walk  # noqa: E402
from exvideo.errors import VideoStageError  # noqa: E402

DEFAULT_WALK = "HYOU1,HSP90B1,CDC37,SRC"


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Build the walk, facts and prompt; ask Claude unless --backend none; save everything."""
    ap = argparse.ArgumentParser(description="walk -> figure 17 facts -> the team's prompt -> 16 bars")
    ap.add_argument("--walk", default=DEFAULT_WALK)
    ap.add_argument("--start", help="random walk from this node instead of --walk")
    ap.add_argument("--steps", type=int, default=3)
    ap.add_argument("--seed", type=int)
    ap.add_argument("--backend", choices=(*lyrics.BACKENDS, "none"), default="cli")
    ap.add_argument("--model", default=lyrics.DEFAULT_MODEL)
    ap.add_argument("--out", default=os.environ.get("HACK_OUT", str(Path.home() / "Desktop/output/hackathon-2026-track1/network")))
    args = ap.parse_args(argv)
    try:
        out = Path(args.out); net = network.load_network(out)
        if args.start:
            path, seed = walk.random_walk(walk.resolve_node(args.start, net), args.steps, net, args.seed)
        else:
            path, seed = tuple(walk.resolve_node(n, net) for n in args.walk.split(",") if n.strip()), None
            walk.check_walk(path, net)
        facts = network.build_facts(path, net, out); prompt = lyrics.build_prompt(path, facts)
        dest = out / "video" / "_".join(path); dest.mkdir(parents=True, exist_ok=True)
        (dest / "walk_facts.json").write_text(json.dumps(facts, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (dest / "lyrics_prompt.md").write_text(prompt + "\n", encoding="utf-8")
        print(f"walk {' -> '.join(path)}" + (f" (seed {seed})" if seed is not None else "") + f": facts and prompt -> {dest}")
        if args.backend == "none":
            return 0
        lyr = lyrics.write_lyrics(prompt, path, args.backend, args.model)
        (dest / "lyrics_raw.txt").write_text(lyr.pop("raw"), encoding="utf-8")
        (dest / "lyrics.json").write_text(json.dumps({"walk": list(path), "seed": seed, **lyr, "model": args.model, "backend": args.backend}, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (dest / "suno_lyrics.txt").write_text(lyrics.suno_text(lyr, path) + "\n", encoding="utf-8")
        (dest / "suno_style.txt").write_text(lyr["suno_style"] + "\n", encoding="utf-8")
        md = [f"# {lyr['title']}", ""] + [line for i, n in enumerate(path) for line in
              ([f"**{n}: {lyr['personas'].get(n, '')}**", ""] + [b["text"] for b in lyr["bars"] if b["node"] == n] + [""])]
        (dest / "lyrics.md").write_text("\n".join(md), encoding="utf-8")
        print("\n".join(md))
        return 0
    except VideoStageError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
