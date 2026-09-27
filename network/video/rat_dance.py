#!/usr/bin/env python3
# =====================================================================================================
# video/rat_dance.py — THE DANCING RAT, ON THE BEAT OF ANY SONG
# =====================================================================================================
#
# PURPOSE
#   Give it any song; it makes a video of the rat-dance meme stepping exactly on the song's beats:
#     1. finds the song's beats (exvideo.beats: onsets -> tempo -> dynamic-programming beat tracking);
#     2. cleans the rat GIF (exvideo.dancer: drops the solid-green flash frame, keeps ONE colour clip, removes the
#        green fringe) and finds its steps (the bottom of each bob) and a seamless loop; dresses it in the red Team
#        2-PAC bandana, knot in front (exvideo.costume; --no-bandana for the plain rat);
#     3. for every video frame, picks the rat frame so that step k's hit lands exactly on beat k (exvideo.dancer.frame_at,
#        the same rule the music-video renderer uses) — the dance speeds up and slows down with the song, never drifts;
#     4. ffmpeg scales the rat up, puts it on a background (or keeps it transparent) and adds the song.
#   The music-video pipeline (make_music_video.py) draws the same rat at the side of its fly-through.
#
# HOW TO RUN
#   python3 network/video/rat_dance.py ~/Downloads/song.mp3                          # -> song_rat.mp4 next to the song
#   python3 network/video/rat_dance.py song.mp3 --size 1920x1080 --bg black --out ~/Desktop/rat.mp4
#   python3 network/video/rat_dance.py song.mp3 --transparent                        # ProRes 4444 .mov, for editing
#   python3 network/video/rat_dance.py song.mp3 --bpm 170                            # if it dances at half/double time
#   python3 network/video/rat_dance.py song.mp3 --nudge -0.03                        # shift the steps earlier by 30 ms
#   Needs: ffmpeg (brew install ffmpeg), numpy, Pillow (bash network/video/setup.sh; switches to .venv by itself); the rat GIF.
#
# OUTPUTS: the video (default <song>_rat.mp4 beside the song) and <video>.beats.json (tempo, beats, dancer loop)
#
# KNOWN LIMITS
#   One tempo is searched for per song (the beat tracker follows drift, not a sudden tempo change). Stretches with no
#   beat (a silent intro) still get a steady beat, continued from the nearest detected beats. If the tracker locks
#   onto half or double time the rat still steps on the beat, just at a different rate; --bpm or --steps-per-beat fix it.
# =====================================================================================================
from __future__ import annotations

import argparse
import json
import os
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional, Sequence

HERE = Path(__file__).resolve().parent
# Use the video tools' own environment (network/video/.venv, made by setup.sh) when this Python lacks numpy or Pillow.
_VENV = HERE / ".venv" / "bin" / "python"
try:
    import numpy as np
    import PIL  # noqa: F401
except ImportError:
    if _VENV.exists() and os.environ.get("EXVIDEO_REEXEC") != "1":
        os.environ["EXVIDEO_REEXEC"] = "1"; os.execv(str(_VENV), [str(_VENV), str(Path(__file__).resolve()), *sys.argv[1:]])
    raise
sys.path.insert(0, str(HERE))
from exvideo import audio, beats, dancer  # noqa: E402
from exvideo.errors import VideoStageError  # noqa: E402

DEFAULT_DANCER = Path(os.environ.get("HACK_EXT_DANCER", str(Path.home() / "Desktop/output/hackathon-2026-track1/external/dancer/rat_dance_transparent.gif")))
BACKGROUNDS = {"white": "white", "black": "black", "green": "0x00FF00"}     # green = chroma-key green


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(description="The dancing rat, stepping on the beat of any song.")
    ap.add_argument("song", help="the song (any format ffmpeg reads: mp3, m4a, wav, ...)")
    ap.add_argument("--out", help="the video (default: <song>_rat.mp4, or .mov with --transparent, beside the song)")
    ap.add_argument("--dancer", default=str(DEFAULT_DANCER), help="the dancer GIF (default: %(default)s)")
    ap.add_argument("--clip", type=int, default=0, help="which single-colour clip of the GIF (the rat: 0 grey, 1 red, 2 teal)")
    ap.add_argument("--no-bandana", action="store_true", help="the plain rat, without the red Team 2-PAC bandana")
    ap.add_argument("--steps-per-beat", type=float, choices=(0.5, 1.0, 2.0), help="rat steps per beat (default: closest to the rat's own speed)")
    ap.add_argument("--bpm", type=float, help="the song's rough tempo, if the rat dances at half or double time (search within x1.25 of it)")
    ap.add_argument("--nudge", type=float, default=0.0, help="shift every step by this many seconds (negative = earlier)")
    ap.add_argument("--size", default="720x1280", help="video size WxH (default %(default)s, portrait)")
    ap.add_argument("--rat-height", type=float, default=0.8, help="rat height as a fraction of the video height (default %(default)s)")
    ap.add_argument("--bg", default="white", help="background: white, black, green (chroma key) or any ffmpeg colour (default %(default)s)")
    ap.add_argument("--transparent", action="store_true", help="no background: ProRes 4444 .mov with alpha, for video editors")
    ap.add_argument("--fps", type=int, default=30)
    args = ap.parse_args(argv)
    try:
        song = audio.clean_path(args.song)
        if not song.exists():
            raise VideoStageError(f"song not found: {song}")
        out = Path(args.out).expanduser() if args.out else song.with_name(song.stem + "_rat" + (".mov" if args.transparent else ".mp4"))
        W, H = (int(x) for x in args.size.lower().split("x"))
        t0 = time.time()

        # 1. the beats
        total = audio.duration(song)
        bt = beats.detect_beats(song, bpm_hint=args.bpm)
        if len(bt["beats"]) < 2:
            raise VideoStageError("no beat found in the song")
        period = statistics.median(b - a for a, b in zip(bt["beats"], bt["beats"][1:]))
        grid = [b + args.nudge for b in beats.extend_grid(bt["beats"], total)]
        print(f"[1/3] beats: {bt['bpm']} BPM, {len(bt['beats'])} beats over {total:.0f} s")

        with tempfile.TemporaryDirectory() as tmp:
            # 2. the rat
            dz = dancer.prepare(Path(args.dancer).expanduser(), Path(tmp) / "frames", clip=args.clip,
                                steps_per_beat=args.steps_per_beat, beat_period=period, bandana=not args.no_bandana)
            loop = dancer.read_rgba_pngs(dz["frames"])
            print(f"[2/3] rat: clip {dz['clip']} (GIF frames {dz['gif_frames'][0]}-{dz['gif_frames'][1]}), {dz['steps_per_loop']} steps per loop, "
                  f"{dz['steps_per_beat']:g} step(s) per beat, {dz['dropped_flash_frames']} flash frame(s) dropped")

            # 3. frame by frame: the rat at native size, piped to ffmpeg, which scales, places and adds the song
            n = int(np.ceil(total * args.fps)); h, w = loop.shape[1:3]
            rat_h = int(round(H * args.rat_height / 2)) * 2
            bg = BACKGROUNDS.get(args.bg, args.bg)
            if args.transparent:
                graph = f"color=c=black@0.0:s={W}x{H}:r={args.fps},format=rgba[bg];[0:v]scale=-2:{rat_h}:flags=lanczos[r];[bg][r]overlay=(W-w)/2:(H-h)/2:shortest=1,format=yuva444p10le[v]"
                codec = ["-c:v", "prores_ks", "-profile:v", "4444", "-pix_fmt", "yuva444p10le"]
            else:
                graph = f"color=c={bg}:s={W}x{H}:r={args.fps}[bg];[0:v]scale=-2:{rat_h}:flags=lanczos[r];[bg][r]overlay=(W-w)/2:(H-h)/2:shortest=1,format=yuv420p[v]"
                codec = ["-c:v", "libx264", "-preset", "medium", "-crf", "18", "-pix_fmt", "yuv420p"]
            out.parent.mkdir(parents=True, exist_ok=True)
            cmd = ["ffmpeg", "-v", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgba", "-s", f"{w}x{h}", "-r", str(args.fps), "-i", "-",
                   "-i", str(song), "-filter_complex", graph, "-map", "[v]", "-map", "1:a", *codec,
                   "-c:a", "pcm_s16le" if args.transparent else "aac", *([] if args.transparent else ["-b:a", "256k"]), "-shortest", str(out)]
            ff = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            for f in range(n):
                ff.stdin.write(loop[dancer.frame_at(f / args.fps, grid, dz["hits"], dz["n_frames"], dz["steps_per_beat"], args.fps)].tobytes())
            ff.stdin.close(); err = ff.stderr.read().decode(); ff.wait()
            if ff.returncode != 0:
                raise VideoStageError(f"ffmpeg failed: {err[:500]}")

        meta = {"song": str(song), "bpm": bt["bpm"], "beats": bt["beats"], "nudge": args.nudge,
                "dancer": {k: v for k, v in dz.items() if k != "frames"}, "fps": args.fps, "size": [W, H]}
        out.with_suffix(".beats.json").write_text(json.dumps(meta, indent=1) + "\n", encoding="utf-8")
        print(f"[3/3] video: {out}  ({n} frames, {time.time() - t0:.0f} s)")
        return 0
    except VideoStageError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
