#!/usr/bin/env python3
# =====================================================================================================
# video/01_lyrics_from_walk.py — MUSIC VIDEO, PART 1: you name a node -> the team's random walk -> 16 bars of lyrics
# =====================================================================================================
#
# PURPOSE
#   The first half of make_music_video.py, for when the song is made later (or by someone else):
#     1. asks for a node / feature of the joint network (a gene symbol or metabolite name, e.g. HYOU1);
#     2. walks 3 steps from it with the TEAM's walker (random_walk/random_walks.R: step probabilities from the chosen
#        arm's edge weights, EE = endurance by default), every step re-checked against the physical edges;
#     3. gathers the figure 17 facts of the 4 nodes and 3 edges, asks Claude for 16 bars with the team's prompt;
#     4. puts the Suno-ready lyrics on the clipboard, prints the Suno style, opens suno.com/create, and STOPS,
#        printing the command for part 2 (02_video_from_song.py), which takes the song back and makes the video.
#
# HOW TO RUN
#   python3 network/video/01_lyrics_from_walk.py                        # asks for the node
#   python3 network/video/01_lyrics_from_walk.py --start HYOU1 --arm RE --seed 7
#   python3 network/video/01_lyrics_from_walk.py --walk HYOU1,HSP90B1,CDC37,SRC      # a fixed walk
#   Any make_music_video.py option works here too (--walker builtin, --model, --no-open, ...).
#
# OUTPUTS ($HACK_OUT/video/<walk>/): walk.json, walk_facts.json, lyrics_prompt.md, lyrics_raw.txt, lyrics.json,
#   lyrics.md, suno_lyrics.txt, suno_style.txt — the folder part 2 reads.
# =====================================================================================================
from __future__ import annotations

import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
_VENV = HERE / ".venv" / "bin" / "python"
try:                                                   # the video tools' environment (setup.sh), as make_music_video.py does
    import faster_whisper  # noqa: F401
except ImportError:
    if _VENV.exists() and os.environ.get("EXVIDEO_REEXEC") != "1":
        os.environ["EXVIDEO_REEXEC"] = "1"; os.execv(str(_VENV), [str(_VENV), str(Path(__file__).resolve()), *sys.argv[1:]])
sys.path.insert(0, str(HERE))
import make_music_video  # noqa: E402

if __name__ == "__main__":
    sys.exit(make_music_video.main(["--stop-after-lyrics", *sys.argv[1:]]))
