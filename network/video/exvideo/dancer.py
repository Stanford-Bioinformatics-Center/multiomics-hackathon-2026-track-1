"""The dancer: a looping GIF whose dance is retimed onto the song's beats.

The GIF is split into frames (ffmpeg, transparency kept). Its own step rhythm is measured from how much the picture
changes frame to frame (the autocorrelation of that motion curve gives the step period), so we know how many steps
one loop holds. In the video, the loop then advances by BEATS, not seconds: each dance step lands on a beat of the
song, whatever its tempo.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict

import numpy as np

from .errors import InputError


def extract_frames(gif: Path, dest: Path) -> int:
    """Write the GIF's frames as PNG files dest/0001.png ... (alpha kept); returns the number of frames."""
    if not gif.exists():
        raise InputError(f"dancer GIF not found: {gif}")
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(gif), "-fps_mode", "passthrough", str(dest / "%04d.png")], capture_output=True, text=True)
    n = len(list(dest.glob("*.png")))
    if res.returncode != 0 or n < 2:
        raise InputError(f"could not read the frames of {gif}: {res.stderr[:300]}")
    return n


def steps_per_loop(gif: Path, n_frames: int) -> int:
    """How many dance steps (= song beats) one loop of the GIF holds.

    The dancer's side-to-side sway is measured from its SILHOUETTE (the alpha mask, so the GIF's colour changes do not
    count): the horizontal centre of mass per frame, smoothed. Its autocorrelation over realistic sway periods (a
    sixth to half of the loop) gives the sway period; one full sway (left and right) is two steps. The result is
    rounded to an even number of steps so the loop wraps cleanly; at least 2.
    """
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(gif), "-fps_mode", "passthrough", "-vf", "format=rgba,alphaextract,scale=48:-1,format=gray",
                          "-f", "rawvideo", "-"], capture_output=True)
    raw = np.frombuffer(res.stdout, dtype=np.uint8)
    if raw.size == 0 or raw.size % n_frames:
        return 8
    f = raw.reshape(n_frames, -1).astype(float); width = 48
    cols = np.arange(width)
    cx = np.array([(fr.reshape(-1, width) * cols).sum() / max(fr.sum(), 1.0) for fr in f])      # horizontal centre of mass
    k = 5; cx = np.convolve(np.r_[cx[-k:], cx, cx[:k]], np.ones(k) / k, mode="same")[k:-k]        # smooth, wrapping
    m = cx - cx.mean()
    ac = np.array([np.dot(m, np.roll(m, lag)) for lag in range(n_frames)])
    lags = np.arange(max(3, n_frames // 12), n_frames // 2 + 1)
    period = float(lags[np.argmax(ac[lags])])
    steps = 2.0 * n_frames / period
    return max(2, int(2 * round(steps / 2)))


def prepare(gif: Path, dest: Path, beats_per_step: int = 1) -> Dict[str, object]:
    """Frames on disk + how many song beats one loop should take."""
    n = extract_frames(gif, dest)
    steps = steps_per_loop(gif, n)
    probe = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries", "stream=width,height", "-of", "csv=p=0", str(gif)], capture_output=True, text=True)
    try:
        w, h = (int(x) for x in probe.stdout.strip().split(",")[:2]); aspect = w / h
    except ValueError:
        aspect = 0.6
    return {"frames": [str(p) for p in sorted(dest.glob("*.png"))], "n_frames": n, "steps_per_loop": steps, "beats_per_loop": steps * beats_per_step, "aspect": round(aspect, 4)}
