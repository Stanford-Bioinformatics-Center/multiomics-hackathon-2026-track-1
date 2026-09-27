"""The lyrics: the team's prompt (verbatim) + the figure 17 facts -> Claude -> 16 checked bars (+ Suno text)."""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any, Callable, Dict, List, Sequence, Tuple

from .errors import ModelError

BARS_PER_NODE = 4
DEFAULT_MODEL = "claude-opus-5-5"

# The team's prompt, verbatim; {walk} is replaced by the walk, e.g. "HYOU1 -> HSP90B1 -> CDC37 -> SRC".
PROMPT_TEMPLATE = (
    "give me 16 bars of 2pac rap lyrics which summarize the most interesting story from these 4 nodes: {walk}. "
    "4 bars per node, and write the lyrics based on what is most interesting based on the figure 17 graph data. "
    "right now, we will keep the walk listed ({walk}) and later we will extend the script to do a random walk on a "
    "prompted node"
)

# The team's positive control: lyrics they loved, for the HYOU1 -> HSP90B1 -> CDC37 -> SRC walk. Used ONLY as a
# style reference (voice, structure, how data become lines); the model is told not to reuse its lines.
STYLE_EXAMPLE = """I. HYOU1: the lookout
Four hours after the whistle, I'm the first one on the block,
Endurance or the iron, either way my numbers pop,
Thirty-somethin' percent up, and the muscle hear me knock,
Only three lines on my phone, but I'm the spark that make it rock.

II. HSP90B1: the stash house (GRP94)
Six sugars on my jacket, N-linked, I stay dressed,
Ten phosphate spots on file, but none of 'em confessed,
Me and HYOU1 tight on the long run, the endurance road the best,
Diabetes checked my papers and I passed the test.

III. CDC37: the middleman
I'm the one who walks the kinases in through the door,
YES1 hit my line, now I'm movin' more and more,
Endurance gave me half, but the iron made me soar,
Link to SRC doubled up when the weights hit the floor.

IV. SRC: the kingpin
Twelve roads lead to me, I'm the hub of the whole map,
Twenty-seven on my payroll: paxillin, integrin, STAT,
Nine-forty-seven sites I can hit, and that's a fact,
Resistance crowned me king, but phosphatases keep me in check."""

PROMPT_SUFFIX = """

Write ORIGINAL lyrics in the style of 1990s West Coast rap (do not quote or reuse lines from any existing song, and do
not reuse the lines of the style example below). Give every node a street persona that fits its biology (e.g. a hub is
"the kingpin", a chaperone that hands proteins on is "the middleman"). Turn the most interesting facts into lines:
strong or arm-specific edges, the hub, disease and ageing directions, phosphosites or sugars, who is on its line. Keep
every biological statement true to the data; round numbers the way a rapper would say them. Bars 1-4 are about node 1,
bars 5-8 node 2, bars 9-12 node 3, bars 13-16 node 4; hand off to the next node where it fits. Each bar is one line.

STYLE EXAMPLE (the voice we want; a different walk may have different facts):
{style}

FIGURE 17 GRAPH DATA (JSON):
{facts}

Answer with JSON only, exactly this shape (no other text):
{{"title": "<song title>",
  "suno_style": "<Suno style prompt, max 120 characters: genre, beat, tempo, voice; NO artist names>",
  "sections": [{{"node": "<node 1>", "persona": "<persona>", "bars": ["<bar 1>", "<bar 2>", "<bar 3>", "<bar 4>"]}}, ... one per node, in walk order]}}
"""


BANNED_FILE = Path(__file__).resolve().parent / "suno_banned.txt"
MAX_WORDS = 12                  # per bar: one rap bar at ~92 BPM is ~2.6 s, so 16 short bars sing in ~45 s and the
                                # whole Suno song (intro, outro, breaks) stays well under 1:30

LENGTH_RULES = """

LENGTH RULES (the song must stay under 1:30):
- Every bar is ONE short line of at most {max_words} words — never two lines joined by a comma.
- The persona is 2-4 words (e.g. "the kingpin"), no explanation in brackets.
- The Suno style must include "short intro"."""

SUNO_RULES = """

SUNO RULES (Suno refuses lyrics that name an artist, producer, label or DJ, or contain a producer tag):
- Name no real musician, producer, label, DJ or crew anywhere (bars, title, personas, style), not even as a nickname.
- Never use these words, in any form or plural (Suno has rejected them or they are artist names): {banned}.
  The style example above says "phosphate" — Suno rejected that word; say "phospho mark", "phospho spot" or "P-tag".
- Prefer plain science words; if a word could double as a stage name or producer tag, pick another."""


def banned_terms() -> List[str]:
    """The Suno-banned terms (exvideo/suno_banned.txt), lower case, comments stripped."""
    out = []
    for line in BANNED_FILE.read_text(encoding="utf-8").splitlines():
        t = line.split("#", 1)[0].strip().lower()
        if t:
            out.append(t)
    return out


def add_banned(term: str) -> bool:
    """Add a word Suno rejected to the list (dated); False if it was already there."""
    t = term.strip().lower()
    if not t or t in banned_terms():
        return False
    from datetime import date
    with BANNED_FILE.open("a", encoding="utf-8") as fh:
        fh.write(f"{t:<23} # rejected by Suno, added {date.today().isoformat()}\n")
    return True


def _banned_rx(terms: Sequence[str]):
    """One regex for all terms: whole words, case-insensitive, plural s / es allowed."""
    alts = sorted((re.escape(t) for t in terms), key=len, reverse=True)
    return re.compile(r"(?<![\w])(" + "|".join(alts) + r")(?:e?s)?(?![\w])", re.I) if alts else None


def find_banned(lyr: Dict[str, Any], terms: Sequence[str] = ()) -> List[Tuple[str, str]]:
    """Every place a banned term appears: (where, term) with where = 'bar 7', 'title', 'persona NODE' or 'style'."""
    rx = _banned_rx(list(terms) or banned_terms())
    if rx is None:
        return []
    hits = []
    fields = [("title", lyr.get("title", "")), ("style", lyr.get("suno_style", ""))]
    fields += [(f"persona {n}", p) for n, p in lyr.get("personas", {}).items()]
    fields += [(f"bar {b['bar']}", b["text"]) for b in lyr.get("bars", [])]
    for where, text in fields:
        hits += [(where, m.group(0)) for m in rx.finditer(text or "")]
    return hits


def find_long(lyr: Dict[str, Any], max_words: int = MAX_WORDS) -> List[Tuple[str, str]]:
    """Every bar longer than `max_words` words, and every persona longer than 5 words: (where, 'N words')."""
    out = [(f"bar {b['bar']}", f"{len(b['text'].split())} words") for b in lyr.get("bars", []) if len(b["text"].split()) > max_words]
    out += [(f"persona {n}", f"{len(p.split())} words") for n, p in lyr.get("personas", {}).items() if len(p.split()) > 5]
    return out


def build_prompt(walk: Sequence[str], facts: Dict[str, Any]) -> str:
    """The full prompt: the team's words, then the originality / style instructions, the example, the data, and the
    Suno rules (with the current banned-word list)."""
    return (PROMPT_TEMPLATE.format(walk=" -> ".join(walk)) + PROMPT_SUFFIX.format(style=STYLE_EXAMPLE, facts=json.dumps(facts, indent=2, ensure_ascii=False))
            + SUNO_RULES.format(banned=", ".join(banned_terms())) + LENGTH_RULES.format(max_words=MAX_WORDS))


def call_cli(prompt: str, model: str) -> str:
    """Ask Claude through the Claude Code command line (non-interactive: claude -p)."""
    exe = shutil.which("claude")
    if exe is None:
        raise ModelError("Claude Code ('claude') is not on the PATH; install it or use --backend api")
    try:
        res = subprocess.run([exe, "-p", prompt, "--model", model, "--output-format", "text"], capture_output=True, text=True, timeout=600)
    except subprocess.TimeoutExpired as err:
        raise ModelError("the Claude Code call timed out (600 s)") from err
    if res.returncode != 0:
        raise ModelError(f"claude -p failed (exit {res.returncode}): {res.stderr.strip()[:500]}")
    return res.stdout


def call_api(prompt: str, model: str) -> str:
    """Ask Claude through the Anthropic Messages API (needs ANTHROPIC_API_KEY)."""
    key = os.environ.get("ANTHROPIC_API_KEY")
    if not key:
        raise ModelError("ANTHROPIC_API_KEY is not set; use --backend cli")
    body = json.dumps({"model": model, "max_tokens": 2500, "messages": [{"role": "user", "content": prompt}]}).encode()
    req = urllib.request.Request("https://api.anthropic.com/v1/messages", data=body, method="POST",
                                 headers={"x-api-key": key, "anthropic-version": "2023-06-01", "content-type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            out = json.loads(r.read().decode())
    except urllib.error.HTTPError as err:
        raise ModelError(f"Anthropic API error {err.code}: {err.read().decode()[:500]}") from err
    except urllib.error.URLError as err:
        raise ModelError(f"could not reach the Anthropic API: {err.reason}") from err
    return "".join(b.get("text", "") for b in out.get("content", []) if b.get("type") == "text")


BACKENDS: Dict[str, Callable[[str, str], str]] = {"cli": call_cli, "api": call_api}


def parse_lyrics(text: str, walk: Sequence[str]) -> Dict[str, Any]:
    """Pull the JSON out of the answer and check it: one section per node in walk order, 4 non-empty bars each."""
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        raise ModelError("the answer contains no JSON object")
    try:
        obj = json.loads(m.group(0))
    except json.JSONDecodeError as err:
        raise ModelError(f"the answer's JSON does not parse: {err}") from err
    sections = obj.get("sections")
    if not isinstance(sections, list) or len(sections) != len(walk):
        raise ModelError(f"expected {len(walk)} sections (one per node), got {len(sections) if isinstance(sections, list) else 'none'}")
    bars: List[Dict[str, Any]] = []; personas: Dict[str, str] = {}
    for i, (sec, node) in enumerate(zip(sections, walk)):
        if not isinstance(sec, dict) or sec.get("node") != node:
            raise ModelError(f"section {i + 1} is about {sec.get('node') if isinstance(sec, dict) else '?'!r} but should be {node!r} (walk order)")
        lines = sec.get("bars")
        if not isinstance(lines, list) or len(lines) != BARS_PER_NODE or not all(isinstance(x, str) and x.strip() for x in lines):
            raise ModelError(f"section {i + 1} ({node}) must have exactly {BARS_PER_NODE} non-empty bars")
        personas[node] = str(sec.get("persona") or "").strip()
        start = len(bars)                                  # bar numbers run 1-16 across the sections
        bars += [{"bar": start + j + 1, "node": node, "text": x.strip()} for j, x in enumerate(lines)]
    style = str(obj.get("suno_style") or "90s west coast hip hop, g-funk, deep bass, laid-back groove, male rap vocal, 92 bpm")[:200]
    return {"title": str(obj.get("title") or "Untitled").strip(), "suno_style": style, "personas": personas, "bars": bars}


def write_lyrics(prompt: str, walk: Sequence[str], backend: str, model: str) -> Dict[str, Any]:
    """Ask the model, check the answer, ask once more (saying what was wrong) if it is malformed."""
    if backend not in BACKENDS:
        raise ModelError(f"unknown backend {backend!r}; use one of {', '.join(BACKENDS)}")
    ask = BACKENDS[backend]
    raw = ask(prompt, model)
    try:
        lyr = {**parse_lyrics(raw, walk), "raw": raw}
    except ModelError as first:
        raw = ask(prompt + f"\n\nYour previous answer was rejected: {first}. Return only the JSON, exactly as specified.", model)
        lyr = {**parse_lyrics(raw, walk), "raw": raw}
    return scrub(lyr, walk, backend, model)


def scrub(lyr: Dict[str, Any], walk: Sequence[str], backend: str, model: str, tries: int = 3) -> Dict[str, Any]:
    """Make the lyrics Suno-ready: while a bar (or title / persona / style) contains a banned word, or a bar is longer
    than MAX_WORDS words (or a persona longer than 5), ask Claude to rewrite ONLY those parts, keeping meaning, facts,
    rhyme and rhythm; every other bar is kept word for word."""
    ask = BACKENDS[backend]
    for _ in range(tries):
        hits, longs = find_banned(lyr), find_long(lyr)
        if not hits and not longs:
            return lyr
        bad = sorted({w for w, _ in hits} | {w for w, _ in longs}, key=lambda w: (not w.startswith("bar"), int(w.split()[1]) if w.startswith("bar ") else 0, w))
        words = sorted({t.lower() for _, t in hits})
        current = {"title": lyr["title"], "suno_style": lyr["suno_style"],
                   "sections": [{"node": n, "persona": lyr["personas"].get(n, ""), "bars": [b["text"] for b in lyr["bars"] if b["node"] == n]} for n in walk]}
        why = []
        if hits:
            why.append(f"they contain the Suno-banned word(s) {', '.join(repr(w) for w in words)} in: {', '.join(sorted({w for w, _ in hits}))} "
                       f"(banned: {', '.join(banned_terms())}; for 'phosphate' say 'phospho mark', 'phospho spot' or 'P-tag')")
        if longs:
            why.append(f"these are too long ({', '.join(f'{w}: {n}' for w, n in longs)}): every bar must be ONE line of at most "
                       f"{MAX_WORDS} words and a persona 2-4 words, so the song stays under 1:30")
        fix = ("These rap lyrics are going to Suno, which rejects artist / producer names and producer tags. "
               + "; and ".join(why).capitalize() + ". Rewrite ONLY those parts, keeping the meaning, the facts, rhyme and rhythm. "
               "Change nothing else. Answer with the full JSON in the same shape, no other text.\n\n" + json.dumps(current, indent=2, ensure_ascii=False))
        raw = ask(fix, model)
        new = parse_lyrics(raw, walk)
        keep = {int(w.split()[1]) for w in bad if w.startswith("bar ")}
        lyr = {**lyr,
               "title": new["title"] if "title" in bad else lyr["title"],
               "suno_style": new["suno_style"] if "style" in bad else lyr["suno_style"],
               "personas": {n: (new["personas"].get(n, p) if f"persona {n}" in bad else p) for n, p in lyr["personas"].items()},
               "bars": [nb if b["bar"] in keep else b for b, nb in zip(lyr["bars"], new["bars"])]}
    hits, longs = find_banned(lyr), find_long(lyr)
    if hits or longs:
        raise ModelError(f"after {tries} rewrites the lyrics still break the Suno rules: banned {hits}, too long {longs}")
    return lyr


def suno_text(lyr: Dict[str, Any], walk: Sequence[str]) -> str:
    """The lyrics formatted for Suno's custom-lyrics box: short section tags, one bar per line."""
    out = ["[Intro]", "(Team 2-PAC)", ""]
    for i, node in enumerate(walk):
        out.append(f"[Verse {i + 1}: {node}]")                    # lean tags: Suno may sing or pause for long ones
        out += [b["text"] for b in lyr["bars"] if b["node"] == node] + [""]
    out += ["[Outro]", "(" + " -> ".join(walk) + ")"]
    return "\n".join(out)
