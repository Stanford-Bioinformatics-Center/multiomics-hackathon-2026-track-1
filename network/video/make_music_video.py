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
import random
import shutil
import statistics
import subprocess
import sys
import time
import webbrowser
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).resolve().parent
# Use the video tools' own environment (network/video/.venv, made by setup.sh) when it exists and this Python lacks
# faster-whisper (needed to sync the lyrics to the song); re-start this script under it once.
_VENV = HERE / ".venv" / "bin" / "python"
try:
    import faster_whisper  # noqa: F401
except ImportError:
    if _VENV.exists() and os.environ.get("EXVIDEO_REEXEC") != "1":
        os.environ["EXVIDEO_REEXEC"] = "1"; os.execv(str(_VENV), [str(_VENV), str(Path(__file__).resolve()), *sys.argv[1:]])
sys.path.insert(0, str(HERE))
from exvideo import align, audio, beats, dancer, lyrics, network, render, starts, walk  # noqa: E402

DEFAULT_DANCER = Path(os.environ.get("HACK_EXT_DANCER", str(Path.home() / "Desktop/output/hackathon-2026-track1/external/dancer/rat_dance_transparent.gif")))
DEFAULT_CREDIT = "Dancing rat: original meme by @ratomilton (TikTok), GIF via Tenor"
from exvideo.errors import VideoStageError, WalkError  # noqa: E402

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]


def ask_start(net: network.Network, out: Path, steps: int = 3) -> str:
    """Offer the story markers (exvideo/starts.py: the top markers of T2D muscle, T2D blood and ageing blood) as a
    numbered menu; take a number, any other figure 17 node that a full walk can start from, or Enter for random."""
    try:
        menu = starts.super_list(out, net)
        starts.write_menu(menu, out / "video" / "start_menu.csv")
    except VideoStageError as err:
        print(f"  (no story menu: {err})"); menu = []
    ok = walk.walkable_starts(net, steps)
    if menu:
        print("\nWhere should the walk start? The strongest markers of our three stories")
        print("(big change in the disease, and the winning exercise arm pushes it back the other way):\n")
        story = None
        for i, e in enumerate(menu, 1):
            first = e["stories"][0]["story"]
            if first != story:
                story = first; print(f"  {story.upper()}")
            print(f"  {i:>2}. {starts.describe(e)}")
    print(f"\nPick a number" + (f" (1-{len(menu)})" if menu else "") + f", or type any of the {len(ok)} figure 17 nodes a {steps + 1}-node walk can start from.")
    while True:
        typed = input("Start node (number, name, or Enter = random): ").strip()
        if not typed:
            node = random.SystemRandom().choice([e["node"] for e in menu] or ok); print(f"  random start: {node}"); return node
        if typed.isdigit() and menu:
            k = int(typed)
            if 1 <= k <= len(menu):
                print(f"  start: {menu[k - 1]['node']}"); return str(menu[k - 1]["node"])
            print(f"  pick 1-{len(menu)}"); continue
        try:
            return walk.resolve_start(typed, net, steps)
        except WalkError as err:
            print(f"  {err}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the six stages; stop with a clear message on any problem."""
    ap = argparse.ArgumentParser(description="Team 2-PAC music video: node -> random walk -> lyrics -> Suno -> video")
    ap.add_argument("--start", help="start node (asked interactively if missing)")
    ap.add_argument("--steps", type=int, default=3, help="random-walk steps (default %(default)s: 4 nodes)")
    ap.add_argument("--seed", type=int, help="random-walk seed (default: random, printed and saved)")
    ap.add_argument("--walk", help="a fixed walk instead, comma-separated (e.g. HYOU1,HSP90B1,CDC37,SRC)")
    ap.add_argument("--walker", choices=("team", "builtin"), default="team",
                    help="team (default): the team's random_walk/random_walks.R; builtin: exvideo.walk.random_walk")
    ap.add_argument("--arm", choices=("coin", "EE", "RE"), default="coin",
                    help="team walker: whose edge weights set the step probabilities: coin (default) = a coin flip picks EE or RE at every step; EE = endurance; RE = resistance")
    ap.add_argument("--backend", choices=tuple(lyrics.BACKENDS), default="cli", help="how to ask Claude (default %(default)s)")
    ap.add_argument("--model", default=lyrics.DEFAULT_MODEL)
    ap.add_argument("--audio", help="the song file (skip waiting)")
    ap.add_argument("--max-seconds", type=float, default=90.0,
                    help="cut the song (and video) at this many seconds, with a 3 s fade-out, if it runs longer (default %(default)s = 1:30; 0 = never)")
    ap.add_argument("--stop-after-lyrics", action="store_true",
                    help="stop once the lyrics are written and on the clipboard (01_lyrics_from_walk.py); make the video later with 02_video_from_song.py")
    ap.add_argument("--lyrics", help="reuse saved lyrics (a lyrics.json from an earlier run or 01_lyrics_from_walk.py) instead of writing new ones; its walk is used")
    ap.add_argument("--downloads", default=str(Path.home() / "Downloads"), help="folder watched for the Suno download")
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--size", default="1280x720", help="video size WxH (default %(default)s)")
    ap.add_argument("--intro", type=float, help="seconds before the first bar (default: 10%% of the song, max 8)")
    ap.add_argument("--outro", type=float, help="seconds after the last bar (default: 10%% of the song, max 8)")
    ap.add_argument("--sync", choices=("whisper", "even"), default="whisper",
                    help="whisper (default): listen to the song and start every line when its first word is sung; even: spread the lines evenly")
    ap.add_argument("--whisper-model", default="small.en", help="faster-whisper model for --sync whisper (default %(default)s; base.en is faster)")
    ap.add_argument("--sung-headers", action="store_true", help="(--sync even only) the song also sings each section header; time a slot for it")
    ap.add_argument("--dancer", default=str(DEFAULT_DANCER), help="a looping dancer GIF (transparent background best) drawn at the side, stepping on the beat (default: %(default)s)")
    ap.add_argument("--no-dancer", action="store_true", help="no dancer")
    ap.add_argument("--no-bandana", action="store_true", help="the plain rat, without the red Team 2-PAC bandana")
    ap.add_argument("--dancer-clip", type=int, default=0, help="which single-colour clip of the GIF to dance (0 = first; the rat GIF: 0 grey, 1 red, 2 teal)")
    ap.add_argument("--dancer-steps-per-beat", type=float, choices=(0.5, 1.0, 2.0), help="dance steps per song beat (default: whichever is closest to the GIF's own speed)")
    ap.add_argument("--dancer-credit", default=DEFAULT_CREDIT, help="credit for the dancer GIF, shown only on the closing card")
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
        saved = None; p_steps = None; step_arms = None
        if args.lyrics:
            saved = json.loads(audio.clean_path(args.lyrics).read_text(encoding="utf-8"))
            if not isinstance(saved.get("bars"), list) or not saved.get("walk"):
                raise VideoStageError(f"{args.lyrics} is not a lyrics.json (needs 'walk' and 'bars')")
            path = tuple(saved["walk"]); walk.check_walk(path, net); seed = saved.get("seed")
            lyrics.parse_lyrics(json.dumps({"sections": [{"node": n, "persona": saved.get("personas", {}).get(n, ""), "bars": [b["text"] for b in saved["bars"] if b["node"] == n]} for n in path]}), path)
        elif args.walk:
            path = tuple(walk.resolve_node(n, net) for n in args.walk.split(",") if n.strip()); walk.check_walk(path, net); seed = None
        else:
            start = walk.resolve_start(args.start, net, args.steps) if args.start else ask_start(net, out, args.steps)
            if args.walker == "team" and args.steps == 3:
                path, seed, p_steps, step_arms = walk.team_walk(start, args.arm, args.seed, out, out / "video" / "_walker")
                walk.check_walk(path, net)                  # every step must be a physical edge (our hard gate)
                if len(path) < 4:
                    raise WalkError(f"{start} sits in a piece of the network too small for a 4-node walk")
            else:
                path, seed = walk.random_walk(start, args.steps, net, args.seed); p_steps = None
        arm_note = (f", team walker, coin flip per step: {', '.join(a for a in step_arms[1:])}" if args.arm == "coin" else f", team walker, {args.arm} weights") if p_steps else ""
        print(f"\n[1/4] walk: {' -> '.join(path)}" + (f"   (seed {seed}{arm_note})" if seed is not None else ""))
        dest = out / "video" / "_".join(path); dest.mkdir(parents=True, exist_ok=True)
        (dest / "walk.json").write_text(json.dumps({"walk": list(path), "seed": seed, "steps": len(path) - 1,
                                                    "walker": ("team (random_walk/random_walks.R)" if p_steps else "builtin") if seed is not None else "fixed",
                                                    "arm": args.arm if p_steps else None, "step_arms": step_arms, "p_step": p_steps}, indent=2) + "\n")

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

        if args.stop_after_lyrics:                          # the two-sitting flow: the song comes back later (02_video_from_song.py)
            rec = {"walk": list(path), "seed": seed, "arm": args.arm if p_steps else None, "step_arms": step_arms, "title": lyr["title"], "suno_style": lyr["suno_style"],
                   "personas": lyr["personas"], "bars": [{"bar": b["bar"], "node": b["node"], "text": b["text"]} for b in lyr["bars"]],
                   "model": args.model, "backend": args.backend}
            (dest / "lyrics.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
            here = os.path.relpath(HERE, Path.cwd())
            print(f"\nLyrics saved in {dest}\nWhen the song is ready, make the video with:\n"
                  f"  python3 {here}/02_video_from_song.py \"{dest}\" <the song file>     (or leave the song out: it asks, or waits for the download)")
            return 0

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
        if args.max_seconds and total > args.max_seconds:              # the 1:30 safety net: fade out and end there
            cut = dest / f"song_first_{int(args.max_seconds)}s.wav"
            res = subprocess.run(["ffmpeg", "-v", "error", "-y", "-i", str(song), "-t", f"{args.max_seconds}",
                                  "-af", f"afade=t=out:st={args.max_seconds - 3:.3f}:d=3", str(cut)], capture_output=True, text=True)
            if res.returncode != 0:
                raise VideoStageError(f"could not cut the song to {args.max_seconds:g} s: {res.stderr[:300]}")
            print(f"      the song is {total:.0f} s; cut to {args.max_seconds:g} s with a 3 s fade-out ({cut.name})")
            song = cut; total = audio.duration(song)
        heads = {n: f"{ROMAN[i]}. {n}: {lyr['personas'].get(n, '')}".rstrip(": ") for i, n in enumerate(path)}
        tl = None
        if args.sync == "whisper":
            try:
                print(f"      syncing the lyrics to the song (Whisper {args.whisper_model}; cached after the first run) ...")
                t = time.time(); words = align.transcribe(song, dest / "transcript.json", args.whisper_model)
                lines = [ln for n in path for ln in ([{"bar": 0, "node": n, "text": heads[n], "header": True}] + [b for b in lyr["bars"] if b["node"] == n])]
                timed, rep = align.sync_lines(lines, words, total)
                (dest / "sync_report.json").write_text(json.dumps(rep, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
                tl = audio.timeline_from_lines(timed, path, total)
                print(f"      synced in {time.time() - t:.0f} s: {rep['matched_words']} of {rep['lyric_words']} lyric words heard; vocals {tl['intro']:.1f}-{total - tl['outro']:.1f} s"
                      + (f"; {len(rep['dropped_headers'])} unsung headers left out" if rep["dropped_headers"] else ""))
            except VideoStageError as err:
                print(f"      could not sync with Whisper ({err}); spreading the lines evenly instead", file=sys.stderr)
        if tl is None:
            tl = audio.timeline(lyr["bars"], path, total, args.intro, args.outro, heads if args.sung_headers else None)
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
        if not args.no_dancer and Path(args.dancer).expanduser().exists():
            t = time.time(); bt = beats.detect_beats(song)
            period = statistics.median(b - a for a, b in zip(bt["beats"], bt["beats"][1:])) if len(bt["beats"]) > 1 else None
            dz = dancer.prepare(Path(args.dancer).expanduser(), dest / "dancer_frames", clip=args.dancer_clip, steps_per_beat=args.dancer_steps_per_beat, beat_period=period, bandana=not args.no_bandana)
            spec["dancer"] = {**dz, "beats": beats.extend_grid(bt["beats"], total), "credit": args.dancer_credit, "aspect": dz["aspect"]}
            (dest / "beats.json").write_text(json.dumps(bt, indent=1) + "\n", encoding="utf-8")
            print(f"      dancer: {bt['bpm']} BPM, {len(bt['beats'])} beats; clip {dz['clip']} of {len(dz['clip_lengths'])} (GIF frames {dz['gif_frames'][0]}-{dz['gif_frames'][1]}), "
                  f"{dz['steps_per_loop']} steps per loop, {dz['steps_per_beat']:g} step(s) per beat, {dz['dropped_flash_frames']} flash frame(s) dropped ({time.time() - t:.0f} s)")
        elif not args.no_dancer:
            print(f"      (no dancer: {args.dancer} not found)")
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
