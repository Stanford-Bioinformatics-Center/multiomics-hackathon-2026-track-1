"""Beat tracking for the dancer: tempo and beat times from the song itself (numpy only).

How it works (plain language): the song is decoded to mono; an "onset strength" curve rises wherever the sound
suddenly gets louder in some frequency band (spectral flux: the positive change of the log energy in 40 log-spaced
bands, frame to frame; log-spaced so the kick's few bass bins count as much as the hi-hat's thousand treble bins);
the tempo is the beat period whose pulse train best lands on those onsets: for every candidate period a comb of
pulses is slid over each ~8 s window of the song (so slow drift is fine) and scores the onset strength ABOVE the song's
average that it hits, per second — a faster level only wins if the extra pulses land on real hits, a slower one loses
the beats it skips. (Plain autocorrelation was dropped: a backbeat or off-beat hi-hats made it pick half time or
1.5-beat lags on ordinary drum grooves.) The search runs over 72-176 BPM, where the tapped beat of most songs sits;
a tempo hint moves that band. The beats themselves are then placed by dynamic programming (Ellis 2007, "Beat tracking by
dynamic programming", J New Music Res 36:51): a chain of onsets that are strong AND about one period apart.
"""
from __future__ import annotations

import subprocess
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .errors import InputError

SR, HOP, NFFT, BANDS = 22050, 512, 2048, 40
ONSET_LAG = 0.019    # s; the flux of a window peaks this long before the attack itself (measured on synthetic drum tracks at 70-174 BPM)


def decode(audio: Path, sr: int = SR) -> np.ndarray:
    """The song as mono float32 samples at `sr` Hz (via ffmpeg)."""
    res = subprocess.run(["ffmpeg", "-v", "error", "-i", str(audio), "-ac", "1", "-ar", str(sr), "-f", "f32le", "-"], capture_output=True)
    if res.returncode != 0 or not res.stdout:
        raise InputError(f"could not decode {audio}: {res.stderr.decode()[:300]}")
    return np.frombuffer(res.stdout, dtype=np.float32)


def onset_strength(y: np.ndarray) -> np.ndarray:
    """Spectral flux per hop (positive log-energy change in 40 log-spaced bands, 30 Hz-11 kHz, averaged over bands so
    the kick counts as much as the hi-hat), with a slow local mean removed, scaled to unit sd."""
    starts = np.arange(0, max(1, len(y) - NFFT), HOP)
    win = np.hanning(NFFT).astype(np.float32)
    frames = np.stack([y[s:s + NFFT] * win for s in starts])
    mag = np.abs(np.fft.rfft(frames, axis=1)) ** 2
    freqs = np.fft.rfftfreq(NFFT, 1.0 / SR); edges = np.geomspace(30.0, SR / 2, BANDS + 1)
    band = np.digitize(freqs, edges) - 1; ok = (band >= 0) & (band < BANDS)
    energy = np.stack([mag[:, ok & (band == b)].sum(axis=1) for b in range(BANDS) if (ok & (band == b)).any()], axis=1)
    logs = np.log1p(1e3 * energy / (energy.mean() + 1e-12))
    flux = np.r_[0.0, np.maximum(0.0, np.diff(logs, axis=0)).mean(axis=1)]
    k = int(SR / HOP)                                      # ~1 s moving average
    local = np.convolve(flux, np.ones(k) / k, mode="same")
    o = np.maximum(0.0, flux - local)
    return o / (o.std() + 1e-9)


def tempo_period(o: np.ndarray, lo_bpm: float = 72, hi_bpm: float = 176, window_s: float = 8.0, step: float = 0.05) -> float:
    """Beat period in hops: the pulse-train period (search step `step` hops) that collects the most above-average onset
    strength per second, each ~`window_s` window at its own best phase."""
    fps = SR / HOP
    om = np.maximum(np.maximum(o, np.r_[0.0, o[:-1]]), np.r_[o[1:], 0.0])     # tolerate +-1 hop around each pulse
    om = om - om.mean()
    win = max(int(window_s * fps), 1); starts = range(0, max(1, len(om) - win // 2), win)
    periods = np.arange(60 * fps / hi_bpm, 60 * fps / lo_bpm + step, step)
    score = np.zeros(len(periods))
    for s0 in starts:
        seg = om[s0:s0 + win]; n = len(seg)
        if n < 2:
            continue
        for pi, p in enumerate(periods):
            k = np.arange(0, n / p)
            ph = np.arange(0, int(np.ceil(p)))[:, None]
            idx = np.round(ph + k[None, :] * p).astype(int)
            val = np.where(idx < n, seg[np.minimum(idx, n - 1)], 0.0).sum(axis=1)
            score[pi] += val.max() / n
    return float(periods[int(np.argmax(score))])


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


def detect_beats(audio: Path, bpm_hint: Optional[float] = None) -> Dict[str, object]:
    """Tempo (BPM) and beat times (seconds) of a song. `bpm_hint` (the song's rough BPM, if the tracker picks half or
    double time) narrows the search to within a factor 1.25 of it."""
    o = onset_strength(decode(audio))
    p = tempo_period(o) if bpm_hint is None else tempo_period(o, lo_bpm=bpm_hint / 1.25, hi_bpm=bpm_hint * 1.25)
    frames = track(o, p)
    times = ((frames * HOP + NFFT / 2) / SR + ONSET_LAG).round(3).tolist()     # a frame's time is its window CENTRE, not its start
    return {"bpm": round(60 * SR / HOP / p, 1), "beats": times}


def extend_grid(beats: List[float], total: float) -> List[float]:
    """Continue the beat grid before the first and after the last detected beat (so the dancer never stops)."""
    if len(beats) < 2:
        return beats
    period = float(np.median(np.diff(beats)))
    before = list(np.arange(beats[0] - period, -period, -period))[::-1]
    after = list(np.arange(beats[-1] + period, total + period, period))
    return [round(float(x), 3) for x in before + beats + after]
