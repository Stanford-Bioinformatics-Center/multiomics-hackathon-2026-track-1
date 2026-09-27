#!/usr/bin/env python3
# =====================================================================================================
# video/fix_lyrics.py — SUNO REJECTED A WORD: ban it and rewrite only the bars that use it
# =====================================================================================================
#
# PURPOSE
#   Suno refuses lyrics that contain an artist / producer name or a producer tag ("Your lyrics contain producer tag
#   phosphate - we don't reference specific artists"). This script:
#     1. adds the rejected word(s) to exvideo/suno_banned.txt, so every future lyric run avoids them too (the list goes
#        into the prompt, and every answer is checked against it);
#     2. asks Claude to rewrite ONLY the bars (or title / persona) that contain a banned word, keeping the meaning,
#        facts, rhyme and rhythm; all other bars stay word for word;
#     3. saves the fixed lyrics (lyrics.json, suno_lyrics.txt, lyrics.md; the old lyrics.json is kept as
#        lyrics_before_fix.json), prints them for Suno and copies them to the clipboard.
#
# HOW TO RUN
#   python3 network/video/fix_lyrics.py <walk folder> phosphate            # the word Suno named
#   python3 network/video/fix_lyrics.py <walk folder> word1 word2          # several
#   python3 network/video/fix_lyrics.py <walk folder>                      # just re-check against the current list
# =====================================================================================================
from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from exvideo import audio, lyrics  # noqa: E402
from exvideo.errors import VideoStageError  # noqa: E402

ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Ban the word(s) Suno rejected and rewrite only the bars that use them.")
    ap.add_argument("folder", help="the walk folder (with lyrics.json), or the lyrics.json itself")
    ap.add_argument("words", nargs="*", help="the word(s) Suno rejected")
    ap.add_argument("--backend", choices=tuple(lyrics.BACKENDS), default="cli")
    ap.add_argument("--model", default=lyrics.DEFAULT_MODEL)
    args = ap.parse_args(argv)
    try:
        where = audio.clean_path(args.folder)
        f = where / "lyrics.json" if where.is_dir() else where
        if not f.exists():
            raise VideoStageError(f"no lyrics.json in {where}")
        dest = f.parent
        for w in args.words:
            print(f"banned from now on: {w}" if lyrics.add_banned(w) else f"already banned: {w}")
        rec = json.loads(f.read_text(encoding="utf-8")); walk = rec["walk"]
        lyr = {"title": rec.get("title", "Untitled"), "suno_style": rec.get("suno_style", ""), "personas": rec.get("personas", {}),
               "bars": [{"bar": b["bar"], "node": b["node"], "text": b["text"]} for b in rec["bars"]]}
        hits = lyrics.find_banned(lyr)
        if not hits:
            print("no banned word in these lyrics; nothing to rewrite")
        else:
            print("rewriting: " + ", ".join(f"{w} ({t})" for w, t in hits) + " ...")
            before = {b["bar"]: b["text"] for b in lyr["bars"]}
            lyr = lyrics.scrub(lyr, walk, args.backend, args.model)
            shutil.copy(f, dest / "lyrics_before_fix.json")
            for b in lyr["bars"]:
                if b["text"] != before[b["bar"]]:
                    print(f"  bar {b['bar']:>2}  was: {before[b['bar']]}\n          now: {b['text']}")
            # the saved record keeps its other fields (seed, arm, model, timing if any); timing of changed bars is kept
            timed = {b["bar"]: b for b in rec["bars"]}
            rec.update(title=lyr["title"], suno_style=lyr["suno_style"], personas=lyr["personas"],
                       bars=[{**timed[b["bar"]], "text": b["text"]} for b in lyr["bars"]])
            f.write_text(json.dumps(rec, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        suno = lyrics.suno_text(lyr, walk)
        (dest / "suno_lyrics.txt").write_text(suno + "\n", encoding="utf-8"); (dest / "suno_style.txt").write_text(lyr["suno_style"] + "\n", encoding="utf-8")
        md = [f"# {lyr['title']}", "", f"*Walk: {' -> '.join(walk)}*", ""]
        for i, n in enumerate(walk):
            md += [f"**{ROMAN[i]}. {n}: {lyr['personas'].get(n, '')}**", ""] + [b["text"] for b in lyr["bars"] if b["node"] == n] + [""]
        (dest / "lyrics.md").write_text("\n".join(md), encoding="utf-8")
        if shutil.which("pbcopy"):
            subprocess.run(["pbcopy"], input=suno, text=True)
        print(f"\nTitle: {lyr['title']}\nStyle: {lyr['suno_style']}\n\n{suno}\n\n(on the clipboard; saved in {dest})")
        return 0
    except VideoStageError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
