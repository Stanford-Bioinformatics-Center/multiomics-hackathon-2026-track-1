"""Beat tracking for the dancer: tempo and beat times from the song itself (numpy only).

How it works (plain language): the song is decoded to mono; an "onset strength" curve rises wherever the sound
suddenly gets louder in some frequency band (spectral flux: the positive change of the log spectrum, frame to frame);
the tempo is the beat period at which that curve best repeats (autocorrelation, with a mild preference for 70-140
BPM, where hip-hop sits); the beats themselves are then placed by dynamic programming (Ellis 2007, "Beat tracking by
dynamic programming", J New Music Res 36:51): a chain of onsets that are strong AND about one period apart.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List

import numpy as np

from .errors import InputError

SR, HOP, NFFT = 22050, 512, 2048


def decode(audio: Path, sr: int = SR) -> np.ndarray:
    """The song as mono float32 samples at `sr` Hz (via ffmpeg)."""
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(audio), "-ac", "1", "-ar", str(sr), "-f", "f32le", "-"], capture_output=True)
    if res.returncode != 0 or not res.stdout:
        raise InputError(f"could not decode {audio}: {res.stderr.decode()[:300]}")
    return np.frombuffer(res.stdout, dtype=np.float32)


def onset_strength(y: np.ndarray) -> np.ndarray:
    """Spectral flux per hop (positive log-spectrum change), with a slow local mean removed, scaled to unit sd."""
    starts = np.arange(0, max(1, len(y) - NFFT), HOP)
    win = np.hanning(NFFT).astype(np.float32)
    frames = np.stack([y[s:s + NFFT] * win for s in starts])
    logs = np.log1p(100.0 * np.abs(np.fft.rfft(frames, axis=1)))
    flux = np.r_[0.0, np.maximum(0.0, np.diff(logs, axis=0)).sum(axis=1)]
    k = int(SR / HOP)                                      # ~1 s moving average
    local = np.convolve(flux, np.ones(k) / k, mode="same")
    o = np.maximum(0.0, flux - local)
    return o / (o.std() + 1e-9)


def tempo_period(o: np.ndarray, lo_bpm: float = 60, hi_bpm: float = 180, centre_bpm: float = 100) -> float:
    """Beat period in hops: the autocorrelation peak, weighted toward `centre_bpm` (log-normal prior)."""
    fps = SR / HOP
    ac = np.correlate(o, o, mode="full")[len(o) - 1:]
    lags = np.arange(int(60 * fps / hi_bpm), int(60 * fps / lo_bpm) + 1)
    bpm = 60 * fps / lags
    w = ac[lags] * np.exp(-0.5 * (np.log2(bpm / centre_bpm) / 0.9) ** 2)
    i = int(np.argmax(w))
    if 0 < i < len(w) - 1:                                 # parabolic refinement of the peak
        a, b, c = w[i - 1], w[i], w[i + 1]; d = 0.5 * (a - c) / (a - 2 * b + c + 1e-12)
        return float(lags[i] + d)
    return float(lags[i])


def track(o: np.ndarray, period: float, tightness: float = 100.0) -> np.ndarray:
    """Beat frames by dynamic programming: maximise onset strength minus a penalty for irregular spacing."""
    n = len(o); score = o.astype(float).copy(); back = -np.ones(n, dtype=int)
    lo, hi = int(round(period / 2)), int(round(2 * period))
    for t in range(lo, n):
        taus = np.arange(max(0, t - hi), t - lo + 1)
        if len(taus) == 0:
            continue
        val = score[taus] - tightness * np.log((t - taus) / period) ** 2
        j = int(np.argmax(val))
        if val[j] > 0:
            score[t] += val[j]; back[t] = taus[j]
    t = int(np.argmax(score[-int(round(period)) - 1:]) + n - int(round(period)) - 1)
    beats = [t]
    while back[t] >= 0:
        t = int(back[t]); beats.append(t)
    return np.array(beats[::-1])


def detect_beats(audio: Path) -> Dict[str, object]:
    """Tempo (BPM) and beat times (seconds) of a song."""
    o = onset_strength(decode(audio))
    p = tempo_period(o)
    frames = track(o, p)
    times = (frames * HOP / SR).round(3).tolist()
    return {"bpm": round(60 * SR / HOP / p, 1), "beats": times}


def extend_grid(beats: List[float], total: float) -> List[float]:
    """Continue the beat grid before the first and after the last detected beat (so the dancer never stops)."""
    if len(beats) < 2:
        return beats
    period = float(np.median(np.diff(beats)))
    before = list(np.arange(beats[0] - period, -period, -period))[::-1]
    after = list(np.arange(beats[-1] + period, total + period, period))
    return [round(float(x), 3) for x in before + beats + after]
