#!/usr/bin/env python3
# =====================================================================================================
# video/make_music_video.py — TEAM 2-PAC MUSIC VIDEO, END TO END: node -> random walk -> lyrics -> Suno -> video
# =====================================================================================================
#
# PURPOSE
#   One interactive command for the whole music video, in about five minutes including the trip to Suno:
#     1. asks for a node / feature of the joint network (a gene symbol or metabolite name, e.g. HYOU1);
#     2. random-walks 3 steps from it along the network's PHYSICAL edges (STRING / Rhea), preferring strongly
#        co-regulated links (probability proportional to max |w_EE|, |w_RE|), never revisiting a node;
#     3. gathers the figure 17 graph facts of the 4 nodes and 3 edges and asks Claude — with the team's prompt,
#        verbatim, and the team's favourite lyrics as the style example — for 16 bars (4 per node) + a Suno style;
#     4. copies the Suno-ready lyrics to the clipboard, prints the style prompt and opens suno.com/create;
#     5. waits for the song (drag the file into the terminal, or it picks up the new download in ~/Downloads);
#     6. renders the video: a fly-through of the figure 17 page along the walk (walked edges gold, arm-specific
#        edges red / blue), each node's persona and fact card, the lyrics bar by bar, the song underneath.
#
# HOW TO RUN
#   python3 network/video/make_music_video.py                      # interactive: asks for the node, waits for the song
#   python3 network/video/make_music_video.py --start HYOU1 --seed 7
#   python3 network/video/make_music_video.py --walk HYOU1,HSP90B1,CDC37,SRC --audio ~/Downloads/song.mp3
#   Needs: the pipeline outputs (bash network/run_all.sh) incl. the figure 17 page; Claude Code (`claude`) or
#   ANTHROPIC_API_KEY; Node.js; ffmpeg (brew install ffmpeg). The renderer installs its Node packages on first use.
#
# OUTPUTS ($HACK_OUT/video/<walk>/): walk.json (walk + seed), walk_facts.json, lyrics_prompt.md, lyrics_raw.txt,
#   lyrics.json (the contract: bars with node and timing), lyrics.md, suno_lyrics.txt, suno_style.txt,
#   render_spec.json, music_video.mp4
#
# KNOWN LIMITS
#   Bars are spread evenly over the vocal part of the song (default: 10% intro / outro, at most 8 s each; adjust with
#   --intro / --outro); Suno's phrasing will not match every bar exactly. The lyrics are a creative summary of the
#   data, not a scientific claim.
# =====================================================================================================
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
import webbrowser
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))
from exvideo import audio, lyrics, network, render, walk  # noqa: E402
from exvideo.errors import VideoStageError, WalkError  # noqa: E402

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]


def ask_start(net: network.Network) -> str:
    """Ask for a start node until the name matches the network."""
    while True:
        try:
            return walk.resolve_node(input("Start node / feature (e.g. HYOU1, SRC, Glutathione): "), net)
        except WalkError as err:
            print(f"  {err}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the six stages; stop with a clear message on any problem."""
    ap = argparse.ArgumentParser(description="Team 2-PAC music video: node -> random walk -> lyrics -> Suno -> video")
    ap.add_argument("--start", help="start node (asked interactively if missing)")
    ap.add_argument("--steps", type=int, default=3, help="random-walk steps (default %(default)s: 4 nodes)")
    ap.add_argument("--seed", type=int, help="random-walk seed (default: random, printed and saved)")
    ap.add_argument("--walk", help="a fixed walk instead, comma-separated (e.g. HYOU1,HSP90B1,CDC37,SRC)")
    ap.add_argument("--backend", choices=tuple(lyrics.BACKENDS), default="cli", help="how to ask Claude (default %(default)s)")
    ap.add_argument("--model", default=lyrics.DEFAULT_MODEL)
    ap.add_argument("--audio", help="the song file (skip waiting)")
    ap.add_argument("--lyrics", help="reuse saved lyrics (a lyrics.json from an earlier run or 01_lyrics_from_walk.py) instead of writing new ones; its walk is used")
    ap.add_argument("--downloads", default=str(Path.home() / "Downloads"), help="folder watched for the Suno download")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--size", default="1280x720", help="video size WxH (default %(default)s)")
    ap.add_argument("--intro", type=float, help="seconds before the first bar (default: 10%% of the song, max 8)")
    ap.add_argument("--outro", type=float, help="seconds after the last bar (default: 10%% of the song, max 8)")
    ap.add_argument("--sung-headers", action="store_true", help="the song also sings each section header (e.g. 'I. HYOU1: the lookout'); time a slot for it")
    ap.add_argument("--no-open", action="store_true", help="do not open Suno / the finished video")
    ap.add_argument("--out", default=os.environ.get("HACK_OUT", str(Path.home() / "Desktop/output/hackathon-2026-track1/network")))
    ap.add_argument("--fig", default=os.environ.get("HACK_FIG", str(Path.home() / "Desktop/output/hackathon")))
    args = ap.parse_args(argv)
    t_start = time.time()
    try:
        out = Path(args.out); net = network.load_network(out)
        page = Path(args.fig) / "17_interactive" / "17a_joint_network.html"
        if not page.exists():
            raise VideoStageError(f"figure 17 page not found: {page} (run: bash network/run_all.sh 17s 17i)")
        width, height = (int(x) for x in args.size.lower().split("x"))

        # 1-2. the walk (from saved lyrics, a fixed walk, or a random walk)
        saved = None
        if args.lyrics:
            saved = json.loads(audio.clean_path(args.lyrics).read_text(encoding="utf-8"))
            if not isinstance(saved.get("bars"), list) or not saved.get("walk"):
                raise VideoStageError(f"{args.lyrics} is not a lyrics.json (needs 'walk' and 'bars')")
            path = tuple(saved["walk"]); walk.check_walk(path, net); seed = saved.get("seed")
            lyrics.parse_lyrics(json.dumps({"sections": [{"node": n, "persona": saved.get("personas", {}).get(n, ""), "bars": [b["text"] for b in saved["bars"] if b["node"] == n]} for n in path]}), path)
        elif args.walk:
            path = tuple(walk.resolve_node(n, net) for n in args.walk.split(",") if n.strip()); walk.check_walk(path, net); seed = None
        else:
            start = walk.resolve_node(args.start, net) if args.start else ask_start(net)
            path, seed = walk.random_walk(start, args.steps, net, args.seed)
        print(f"\n[1/4] walk: {' -> '.join(path)}" + (f"   (seed {seed})" if seed is not None else ""))
        dest = out / "video" / "_".join(path); dest.mkdir(parents=True, exist_ok=True)
        (dest / "walk.json").write_text(json.dumps({"walk": list(path), "seed": seed, "steps": len(path) - 1}, indent=2) + "\n")

        # 3. facts -> prompt -> lyrics
        facts = network.build_facts(path, net, out)
        prompt = lyrics.build_prompt(path, facts)
        (dest / "walk_facts.json").write_text(json.dumps(facts, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        (dest / "lyrics_prompt.md").write_text(prompt + "\n", encoding="utf-8")
        if saved is not None:
            print(f"[2/4] reusing the saved lyrics: {args.lyrics}"); t = time.time()
            lyr = {"title": saved.get("title", "Untitled"), "suno_style": saved.get("suno_style", ""), "personas": saved.get("personas", {}),
                   "bars": [{"bar": b["bar"], "node": b["node"], "text": b["text"]} for b in saved["bars"]]}
        else:
            print(f"[2/4] writing the lyrics with Claude ({args.backend}, {args.model}) ...")
            t = time.time(); lyr = lyrics.write_lyrics(prompt, path, args.backend, args.model)
            (dest / "lyrics_raw.txt").write_text(lyr.pop("raw"), encoding="utf-8")
        suno = lyrics.suno_text(lyr, path)
        (dest / "suno_lyrics.txt").write_text(suno + "\n", encoding="utf-8"); (dest / "suno_style.txt").write_text(lyr["suno_style"] + "\n", encoding="utf-8")
        md = [f"# {lyr['title']}", "", f"*Walk: {' -> '.join(path)} · {args.model}*", ""]
        for i, n in enumerate(path):
            md += [f"**{ROMAN[i]}. {n}: {lyr['personas'].get(n, '')}**", ""] + [b["text"] for b in lyr["bars"] if b["node"] == n] + [""]
        (dest / "lyrics.md").write_text("\n".join(md), encoding="utf-8")
        print(f"      done in {time.time() - t:.0f} s\n\n" + "\n".join(md[4:]))

        # 4. Suno (skipped when the song already exists: saved lyrics + a given audio file)
        if saved is not None and args.audio:
            print("[3/4] Suno: skipped (saved lyrics and the song were given)")
        elif shutil.which("pbcopy"):
            subprocess.run(["pbcopy"], input=suno, text=True)
        if not (saved is not None and args.audio):
            print("[3/4] SUNO: the lyrics are on your clipboard (also in suno_lyrics.txt). In Suno, Create -> Custom:")
            print(f"      paste the lyrics; style: {lyr['suno_style']}")
            print(f"      title: {lyr['title']}")
            if not args.no_open:
                webbrowser.open("https://suno.com/create")

        # 5. the song
        if args.audio:
            song = audio.clean_path(args.audio)
        else:
            since = time.time()
            typed = input(f"\n      Drag the song file here and press Enter (or just press Enter to wait for it in {args.downloads}): ")
            if typed.strip():
                song = audio.clean_path(typed)
            else:
                print("      waiting for the download ...")
                song = audio.wait_for_download(Path(args.downloads).expanduser(), since)
        total = audio.duration(song)
        heads = {n: f"{ROMAN[i]}. {n}: {lyr['personas'].get(n, '')}".rstrip(": ") for i, n in enumerate(path)} if args.sung_headers else None
        tl = audio.timeline(lyr["bars"], path, total, args.intro, args.outro, heads)
        rec = {"walk": list(path), "seed": seed, "title": lyr["title"], "suno_style": lyr["suno_style"], "personas": lyr["personas"],
               "audio": str(song), **tl, "model": args.model, "backend": args.backend, "created": datetime.now(timezone.utc).isoformat(timespec="seconds")}
        (dest / "lyrics.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[4/4] rendering the video: {total:.0f} s of song, {args.fps} fps, {width}x{height} ...")

        # 6. the video
        byn = {f["node"]: f for f in facts["nodes"]}
        spec = {"html": str(page), "width": width, "height": height, "fps": args.fps, "duration": total, "title": lyr["title"],
                "subtitle": "a walk through the MoTrPAC exercise network", "walk": list(path), "bars": tl["bars"],
                "segments": [{**s, "label": f"{ROMAN[i]}. {s['node']}", "persona": lyr["personas"].get(s["node"], ""), "facts": network.fact_card(byn[s["node"]])}
                             for i, s in enumerate(tl["segments"])]}
        video = render.render(spec, song, dest)
        print(f"\ndone in {time.time() - t_start:.0f} s: {video}")
        if not args.no_open and shutil.which("open"):
            subprocess.run(["open", str(video)])
        return 0
    except (VideoStageError, KeyboardInterrupt) as err:
        print(f"\nstopped: {err or 'interrupted'}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
