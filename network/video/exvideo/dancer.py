"""The dancer: a looping GIF whose dance steps are pinned onto the song's beats.

What goes wrong with a meme GIF played as-is (plain language):
  - the rat-dance GIF is really several clips glued together, each in a different colour (grey, red, teal), with a
    single solid-green frame between two of them; played straight through, the rat flashes and changes colour;
  - the clip seams are not the same pose, so every loop the rat jumps and seems to stop and restart;
  - the chroma-key cut-out leaves a green fringe around the feet.
What we do instead:
  1. drop "flash" frames (a frame whose opaque area is far larger than usual is a solid-colour card, not the rat);
  2. split the rest into clips wherever the rat's average colour jumps, and keep ONE clip (the first by default);
  3. find the dance's hits: the rat bobs, and the bottom of each bob (the silhouette's lowest point of the centre of
     mass, i.e. the knee bend) is where a step lands, so that is the frame that must sit exactly on a beat;
  4. cut a seamless loop between two hits whose silhouettes match best (an even number of steps, so left/right
     alternation survives the wrap);
  5. remove the green fringe (green pulled down to the larger of red and blue on edge pixels).
The renderer then plays the loop hit-to-hit: hit k is shown exactly on beat k, and the frames between two hits are
spread over the time between the two beats, so the dance speeds up or slows down with the song but never drifts.
"""
from __future__ import annotations

import shutil
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from .errors import InputError


def read_rgba(gif: Path) -> np.ndarray:
    """All frames of the GIF as an array (frames, height, width, 4), uint8, transparency kept."""
    if not gif.exists():
        raise InputError(f"dancer GIF not found: {gif}")
    probe = subprocess.run(["ffprobe", "-v", "error", "-select_streams", "v:0", "-show_entries", "stream=width,height", "-of", "csv=p=0", str(gif)],
                           capture_output=True, text=True)
    try:
        w, h = (int(x) for x in probe.stdout.strip().split(",")[:2])
    except ValueError:
        raise InputError(f"could not read the size of {gif}: {probe.stderr[:300]}")
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(gif), "-fps_mode", "passthrough", "-vf", "format=rgba", "-f", "rawvideo", "-"], capture_output=True)
    raw = np.frombuffer(res.stdout, dtype=np.uint8)
    if res.returncode != 0 or raw.size == 0 or raw.size % (w * h * 4):
        raise InputError(f"could not read the frames of {gif}: {res.stderr.decode()[:300]}")
    f = raw.reshape(-1, h, w, 4)
    if len(f) < 4:
        raise InputError(f"{gif} has too few frames to dance ({len(f)})")
    return f


def clips(rgba: np.ndarray, colour_jump: float = 25.0, flash_area: float = 2.0) -> List[List[int]]:
    """Frame indices grouped into single-colour clips; flash frames (opaque area > flash_area x median) are dropped."""
    alpha = rgba[..., 3] > 128
    area = alpha.sum(axis=(1, 2))
    keep = [i for i in range(len(rgba)) if 0 < area[i] <= flash_area * np.median(area)]
    colour = {i: rgba[i][alpha[i]][:, :3].astype(float).mean(axis=0) for i in keep}
    out: List[List[int]] = []
    for i in keep:
        if out and np.abs(colour[i] - colour[out[-1][-1]]).max() <= colour_jump and i == out[-1][-1] + 1:
            out[-1].append(i)
        else:
            out.append([i])
    return out


def hits(alpha: np.ndarray) -> List[int]:
    """Frames where a step lands: local maxima (lowest points) of the silhouette's vertical centre of mass, smoothed,
    at least half a bob apart."""
    a = alpha.astype(float); rows = np.arange(a.shape[1])
    cy = (a.sum(axis=2) * rows).sum(axis=1) / np.maximum(a.sum(axis=(1, 2)), 1.0)
    cy = np.convolve(np.r_[cy[:1], cy, cy[-1:]], np.ones(3) / 3, mode="same")[1:-1]
    m = cy - cy.mean(); n = len(m)
    ac = np.array([np.dot(m[:n - lag], m[lag:]) for lag in range(n)])
    lo = 4; hi = max(lo + 1, n // 2)
    bob = int(np.argmax(ac[lo:hi]) + lo) if hi > lo else n
    peaks = [i for i in range(1, n - 1) if cy[i] >= cy[i - 1] and cy[i] > cy[i + 1]]
    out: List[int] = []
    for p in peaks:                                        # keep the deeper of two peaks closer than half a bob
        if out and p - out[-1] < bob / 2:
            if cy[p] > cy[out[-1]]:
                out[-1] = p
        else:
            out.append(p)
    return out


def best_loop(alpha: np.ndarray, hit: List[int], tol: float = 0.015) -> Tuple[int, int]:
    """The (start, end) hit indices of the seamless loop: an even number of steps whose end silhouette matches the start.
    Among loops whose mismatch is within `tol` of the best, the one with the most steps wins (less repetitive)."""
    cand = []
    for i in range(len(hit)):
        for j in range(i + 2, len(hit), 2):
            d = float(np.abs(alpha[hit[i]].astype(float) - alpha[hit[j]].astype(float)).mean())
            cand.append((d, j - i, i, j))
    if not cand:
        raise InputError("the dancer clip is too short for a two-step loop")
    best = min(c[0] for c in cand)
    d, _, i, j = max((c for c in cand if c[0] <= best + tol), key=lambda c: (c[1], -c[0]))
    return i, j


def despill(rgba: np.ndarray) -> np.ndarray:
    """Remove the green-screen fringe: on opaque pixels whose green beats both red and blue, cap green at max(red, blue)."""
    out = rgba.copy(); r, g, b = (out[..., k].astype(int) for k in range(3))
    cap = np.maximum(r, b); spill = (g > cap + 8) & (out[..., 3] > 0)
    out[..., 1] = np.where(spill, cap, g).astype(np.uint8)
    return out


def write_pngs(rgba: np.ndarray, dest: Path) -> List[str]:
    """Frames to dest/0001.png ... (via ffmpeg, alpha kept)."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    n, h, w, _ = rgba.shape
    res = subprocess.run(["ffmpeg", "-v", "error", "-f", "rawvideo", "-pix_fmt", "rgba", "-s", f"{w}x{h}", "-i", "-", str(dest / "%04d.png")],
                         input=np.ascontiguousarray(rgba).tobytes(), capture_output=True)
    files = sorted(str(p) for p in dest.glob("*.png"))
    if res.returncode != 0 or len(files) != n:
        raise InputError(f"could not write the dancer frames: {res.stderr.decode()[:300]}")
    return files


def prepare(gif: Path, dest: Path, clip: int = 0, steps_per_beat: Optional[float] = None, beat_period: Optional[float] = None,
            gif_frame_s: float = 0.03) -> Dict[str, object]:
    """Loop frames on disk + where the hits are.

    Returns frames (the loop, in order), hits (frame positions of the steps within the loop, starting at 0; the loop
    wraps back to 0 after the last frame), steps_per_beat (1 = a step on every beat; chosen from 1/2, 1, 2 as the
    rate closest to the GIF's own speed when not given), aspect, and a record of what was dropped.
    """
    rgba = read_rgba(gif)
    groups = clips(rgba)
    if not 0 <= clip < len(groups):
        raise InputError(f"dancer clip {clip} does not exist (the GIF has {len(groups)}: lengths {[len(g) for g in groups]})")
    idx = groups[clip]
    frames = rgba[idx]
    alpha = frames[..., 3] > 128
    hit = hits(alpha)
    if len(hit) < 3:
        raise InputError(f"could not find the dance's steps in clip {clip} of {gif} (found {len(hit)} hits)")
    i, j = best_loop(alpha, hit)
    loop = despill(frames[hit[i]:hit[j]])                  # the end hit is the start hit again, so it is left out
    rel = [h - hit[i] for h in hit[i:j]]
    step_frames = float(np.mean(np.diff(hit[i:j + 1])))
    if steps_per_beat is None:
        natural = step_frames * gif_frame_s                # seconds per step at the GIF's own speed
        steps_per_beat = min((0.5, 1.0, 2.0), key=lambda s: abs(np.log((beat_period or natural) / s / natural)))
    files = write_pngs(loop, dest)
    h, w = loop.shape[1:3]
    return {"frames": files, "n_frames": len(files), "hits": rel, "steps_per_loop": len(rel), "steps_per_beat": steps_per_beat,
            "aspect": round(w / h, 4), "clip": clip, "clip_lengths": [len(g) for g in groups],
            "gif_frames": [int(idx[hit[i]]), int(idx[hit[j]])], "dropped_flash_frames": len(rgba) - sum(len(g) for g in groups)}
