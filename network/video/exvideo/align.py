"""Lyrics <-> song synchronisation: Whisper word timestamps, then a global alignment of sung words to the lyrics.

How it works (plain language):
  1. Whisper (faster-whisper, run locally) transcribes the song with a start / end time for every word it hears.
  2. The known lyrics and the heard words are compared word by word with a global alignment (Needleman-Wunsch):
     two words score by how similar their spelling is after normalising (lower case, punctuation removed, number
     words turned into digits, so "Four" matches "4" and "HYOU1" still partly matches "HU1"); skipping a word
     costs a small penalty. The best alignment tolerates mishearings, extra words and missing words.
  3. Each lyric line starts at the time of its first aligned word (lines with no aligned word are placed between
     their neighbours) and ends when the next line starts. A section header that is not sung (almost none of its
     words heard) is dropped automatically.
The transcript is cached next to the video (keyed by the audio file's size and date), so re-renders are instant.
"""
from __future__ import annotations

import json
import re
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

from .errors import InputError

NUMBER_WORDS = {"zero": "0", "one": "1", "two": "2", "three": "3", "four": "4", "five": "5", "six": "6", "seven": "7", "eight": "8",
                "nine": "9", "ten": "10", "eleven": "11", "twelve": "12", "thirteen": "13", "fourteen": "14", "fifteen": "15",
                "sixteen": "16", "seventeen": "17", "eighteen": "18", "nineteen": "19", "twenty": "20", "thirty": "30", "forty": "40",
                "fifty": "50", "sixty": "60", "seventy": "70", "eighty": "80", "ninety": "90", "hundred": "100", "thousand": "1000"}
GAP = -0.45          # cost of skipping a word on either side
MIN_HEADER_HEARD = 0.34   # a header counts as sung if at least this share of its words is aligned


def normalise(text: str) -> List[str]:
    """Words of a line, lower case, punctuation removed, number words as digits ("thirty-somethin'" -> 30, somethin)."""
    toks = re.findall(r"[a-z0-9]+", text.lower().replace("'", ""))
    return [NUMBER_WORDS.get(t, t) for t in toks]


def similarity(a: str, b: str) -> float:
    """Spelling similarity of two normalised words, 0..1."""
    return 1.0 if a == b else SequenceMatcher(None, a, b).ratio()


def align_words(lyric: Sequence[str], heard: Sequence[str]) -> List[Tuple[int, int]]:
    """Global alignment of lyric words to heard words; returns the matched (lyric index, heard index) pairs."""
    n, m = len(lyric), len(heard)
    score = [[0.0] * (m + 1) for _ in range(n + 1)]; back = [[0] * (m + 1) for _ in range(n + 1)]   # 0 diag, 1 up, 2 left
    for i in range(1, n + 1):
        score[i][0] = i * GAP; back[i][0] = 1
    for j in range(1, m + 1):
        score[0][j] = j * GAP; back[0][j] = 2
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            d = score[i - 1][j - 1] + (2.0 * similarity(lyric[i - 1], heard[j - 1]) - 1.0)
            u = score[i - 1][j] + GAP; l = score[i][j - 1] + GAP
            score[i][j], back[i][j] = max((d, 0), (u, 1), (l, 2))
    pairs: List[Tuple[int, int]] = []; i, j = n, m
    while i > 0 and j > 0:
        if back[i][j] == 0:
            if similarity(lyric[i - 1], heard[j - 1]) >= 0.5:
                pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif back[i][j] == 1:
            i -= 1
        else:
            j -= 1
    return pairs[::-1]


def transcribe(audio: Path, cache: Path, model: str = "small.en") -> List[Dict[str, Any]]:
    """Every heard word with its start / end time (seconds), cached in `cache`."""
    st = audio.stat(); key = {"audio": str(audio), "size": st.st_size, "mtime": int(st.st_mtime), "model": model}
    if cache.exists():
        c = json.loads(cache.read_text(encoding="utf-8"))
        if c.get("key") == key:
            return c["words"]
    try:
        from faster_whisper import WhisperModel
    except ImportError as err:
        raise InputError("faster-whisper is not installed; run: bash network/video/setup.sh (or pip install faster-whisper)") from err
    wm = WhisperModel(model, device="cpu", compute_type="int8")
    segments, _ = wm.transcribe(str(audio), word_timestamps=True, vad_filter=False, beam_size=5)
    words = [{"word": w.word.strip(), "start": round(float(w.start), 3), "end": round(float(w.end), 3)} for s in segments for w in (s.words or [])]
    if not words:
        raise InputError(f"no words were heard in {audio} (is it the right file?)")
    cache.write_text(json.dumps({"key": key, "words": words}, indent=1), encoding="utf-8")
    return words


def sync_lines(lines: Sequence[Dict[str, Any]], words: Sequence[Dict[str, Any]], total: float) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Give every lyric line a start / end from the heard words; drop headers that are not sung.

    `lines`: dicts with "text" (and optionally "header": True), in order. Returns (timed lines, report).
    """
    lyric: List[str] = []; owner: List[int] = []
    for k, ln in enumerate(lines):
        toks = normalise(ln["text"]); lyric += toks; owner += [k] * len(toks)
    heard: List[str] = []; hw: List[int] = []
    for j, w in enumerate(words):
        toks = normalise(w["word"]); heard += toks; hw += [j] * len(toks)
    pairs = align_words(lyric, heard)
    first: Dict[int, float] = {}; last: Dict[int, float] = {}; hits: Dict[int, int] = {}
    for li, hi in pairs:
        k = owner[li]; w = words[hw[hi]]
        first[k] = min(first.get(k, w["start"]), w["start"]); last[k] = max(last.get(k, w["end"]), w["end"]); hits[k] = hits.get(k, 0) + 1
    keep = [k for k, ln in enumerate(lines) if not (ln.get("header") and hits.get(k, 0) < MIN_HEADER_HEARD * max(1, len(normalise(ln["text"]))))]
    starts: Dict[int, Optional[float]] = {k: first.get(k) for k in keep}
    # place lines without any aligned word between their neighbours
    known = [k for k in keep if starts[k] is not None]
    if not known:
        raise InputError("none of the lyrics could be matched to the sung words (different song or lyrics?)")
    for idx, k in enumerate(keep):
        if starts[k] is None:
            prev = next((starts[p] for p in reversed(keep[:idx]) if starts[p] is not None), None)
            nxt = next(((q, starts[q]) for q in keep[idx + 1:] if starts[q] is not None), None)
            if prev is not None and nxt is not None:
                gap_n = keep.index(nxt[0]) - idx + 1; starts[k] = prev + (nxt[1] - prev) * (1 / gap_n)
            elif prev is not None:
                starts[k] = prev + 2.0
            else:
                starts[k] = max(0.0, nxt[1] - 2.0)
    timed: List[Dict[str, Any]] = []
    for idx, k in enumerate(keep):
        end = starts[keep[idx + 1]] if idx + 1 < len(keep) else min(total, (last.get(k) or starts[k] + 2.0) + 0.6)
        timed.append({**lines[k], "start": round(float(starts[k]), 3), "end": round(float(max(end, starts[k] + 0.3)), 3)})
    report = {"lines": len(lines), "kept": len(keep), "dropped_headers": [lines[k]["text"] for k in range(len(lines)) if k not in keep],
              "matched_words": len(pairs), "lyric_words": len(lyric), "heard_words": len(heard),
              "lines_without_match": [lines[k]["text"] for k in keep if k not in first]}
    return timed, report
