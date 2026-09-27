"""The song: wait for the Suno audio, measure it, and lay the 16 bars on its timeline."""
from __future__ import annotations

import shlex
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

from .errors import InputError

AUDIO_EXT = {".mp3", ".wav", ".m4a", ".aac", ".flac", ".ogg"}
PARTIAL_EXT = {".crdownload", ".download", ".part", ".tmp"}


def clean_path(text: str) -> Path:
    """A path typed or dragged into the terminal (quotes and backslash-escaped spaces removed)."""
    text = text.strip()
    parts = shlex.split(text) if text else []
    return Path(parts[0] if parts else text).expanduser()


def wait_for_download(folder: Path, since: float, timeout: float = 900.0, poll: float = 1.0) -> Path:
    """Wait for a new, fully written audio file in `folder` (modified after `since`, size stable for 2 s)."""
    deadline = time.time() + timeout; last: Dict[Path, int] = {}
    while time.time() < deadline:
        if not any(p.suffix.lower() in PARTIAL_EXT for p in folder.iterdir()):
            fresh = sorted((p for p in folder.iterdir() if p.suffix.lower() in AUDIO_EXT and p.stat().st_mtime > since), key=lambda p: p.stat().st_mtime)
            for p in reversed(fresh):
                size = p.stat().st_size
                if size > 0 and last.get(p) == size:
                    return p
                last[p] = size
        time.sleep(poll)
    raise InputError(f"no new audio file appeared in {folder} within {int(timeout)} s")


def duration(audio: Path) -> float:
    """The audio's length in seconds (ffprobe)."""
    if shutil.which("ffprobe") is None:
        raise InputError("ffprobe not found; install ffmpeg (brew install ffmpeg)")
    if not audio.exists():
        raise InputError(f"audio file not found: {audio}")
    res = subprocess.run(["ffprobe", "-v", "error", "-show_entries", "format=duration", "-of", "default=nw=1:nk=1", str(audio)], capture_output=True, text=True)
    try:
        d = float(res.stdout.strip())
    except ValueError as err:
        raise InputError(f"could not read the length of {audio}: {res.stderr.strip()[:300]}") from err
    if d <= 5:
        raise InputError(f"{audio} is only {d:.1f} s long")
    return d


def timeline(bars: Sequence[Dict[str, Any]], walk: Sequence[str], total: float, intro: Optional[float] = None, outro: Optional[float] = None,
             headers: Optional[Dict[str, str]] = None) -> Dict[str, Any]:
    """Spread the bars evenly over the vocal part of the song: [intro, total - outro].

    Default intro / outro = 10% of the song, at most 8 s each (Suno songs open and close with a few bars of beat).
    Every node gets an equal segment; every bar an equal slot within it. With `headers` (node -> header text), the
    song also sings each section header, so every node's segment starts with one extra slot for its header
    (returned among the bars with "header": true).
    """
    intro = min(8.0, 0.10 * total) if intro is None else intro
    outro = min(8.0, 0.10 * total) if outro is None else outro
    if intro + outro >= total - 5:
        raise InputError(f"intro ({intro} s) + outro ({outro} s) leave no time for the verses in a {total:.1f} s song")
    lines: List[Dict[str, Any]] = []
    for n in walk:
        if headers is not None:
            lines.append({"bar": 0, "node": n, "text": headers[n], "header": True})
        lines += [b for b in bars if b["node"] == n]
    per_bar = (total - intro - outro) / len(lines)
    timed: List[Dict[str, Any]] = [{**b, "start": round(intro + i * per_bar, 3), "end": round(intro + (i + 1) * per_bar, 3)} for i, b in enumerate(lines)]
    segs = [{"node": n, "start": min(b["start"] for b in timed if b["node"] == n), "end": max(b["end"] for b in timed if b["node"] == n)} for n in walk]
    return {"duration": round(total, 3), "intro": round(intro, 3), "outro": round(outro, 3), "bars": timed, "segments": segs}
