#!/usr/bin/env python3
# =====================================================================================================
# video/02_video_from_song.py — MUSIC VIDEO, PART 2: the song comes back -> the video
# =====================================================================================================
#
# PURPOSE
#   The second half of make_music_video.py: takes the walk folder that 01_lyrics_from_walk.py wrote and the song made
#   from its lyrics, and renders the video (the lyrics synced to the vocals with Whisper, the fly-through of the
#   figure 17 network along the walk, the dancing rat in its bandana on the beat, the song underneath).
#
# HOW TO RUN
#   python3 network/video/02_video_from_song.py <walk folder> ~/Downloads/song.mp3
#   python3 network/video/02_video_from_song.py <walk folder>            # asks: drag the song in, or Enter to wait
#                                                                        # for the next download to land in ~/Downloads
#   <walk folder> is the folder part 1 printed ($HACK_OUT/video/<walk>/), or its lyrics.json. Any
#   make_music_video.py option works after these (--sync even, --no-dancer, --size 1920x1080, ...).
#
# OUTPUTS (in the walk folder): transcript.json, sync_report.json, beats.json, render_spec.json, music_video.mp4
# =====================================================================================================
from __future__ import annotations

import os
import sys
import time
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
from exvideo import audio  # noqa: E402
from exvideo.errors import VideoStageError  # noqa: E402


def main(argv: list) -> int:
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__ or "", "usage: 02_video_from_song.py <walk folder | lyrics.json> [song] [make_music_video.py options]", sep="")
        return 0 if argv else 2
    try:
        where = audio.clean_path(argv[0])
        lyr = where / "lyrics.json" if where.is_dir() else where
        if not lyr.exists():
            raise VideoStageError(f"no lyrics.json in {where} (run 01_lyrics_from_walk.py first)")
        rest = argv[1:]
        if rest and not rest[0].startswith("-"):
            song, rest = audio.clean_path(rest[0]), rest[1:]
        else:
            since = time.time()
            typed = input("Drag the song file here and press Enter (or just Enter to wait for it in ~/Downloads): ")
            if typed.strip():
                song = audio.clean_path(typed)
            else:
                print("waiting for the download ...")
                song = audio.wait_for_download(Path.home() / "Downloads", since)
        if not song.exists():
            raise VideoStageError(f"song not found: {song}")
    except VideoStageError as err:
        print(f"error: {err}", file=sys.stderr)
        return 2
    return make_music_video.main(["--lyrics", str(lyr), "--audio", str(song), *rest])


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
