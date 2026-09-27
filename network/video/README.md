# Music video (Team 2-PAC)

A walk through the joint network, turned into a rap and a video that flies through the figure 17 page.

| Stage | Script | Does |
|---|---|---|
| 1 | `01_lyrics_from_walk.py` | walk → figure 17 graph facts → the team's prompt → Claude → 16 bars (4 per node) → `lyrics.json` |
| 2 | *(teammate's walk / video script)* | reads `lyrics.json`, flies the figure 17 page node by node with the lyrics on screen, adds the song |

## Stage 1

```bash
python3 network/video/01_lyrics_from_walk.py                  # default walk HYOU1 -> HSP90B1 -> CDC37 -> SRC, via Claude Code (claude -p)
python3 network/video/01_lyrics_from_walk.py --backend none   # only write the facts and the prompt
python3 network/video/01_lyrics_from_walk.py --backend api    # Anthropic API instead (needs ANTHROPIC_API_KEY)
python3 network/video/01_lyrics_from_walk.py --walk A,B,C,D   # any walk along joint-network edges
```

- **The prompt** is the team's, verbatim (`PROMPT_TEMPLATE`), followed by the figure 17 graph data as JSON (per node:
  degree, strength per arm and hub rank, module and pathway, strongest significant exercise response per arm, T2D /
  ageing directions from Öhman 2021 and UK Biobank, responding phosphosites, glycosylation; per edge: STRING score,
  weight per arm, difference, arm specificity) and the output format. The facts and the exact prompt are saved.
- **Checks:** every step of the walk must be a physical edge of the joint network (the network's hard gate); the answer
  must be 16 bars, 4 per node, in walk order (one automatic retry if not). Errors say what is wrong.
- **Random walk (planned):** `--start NODE --steps N --seed S` is reserved; `random_walk()` documents the intended rule
  (step to a neighbour with probability proportional to |w|, no revisits, seeded).

## The contract: `$HACK_OUT/video/<walk>/lyrics.json`

```json
{
  "walk": ["HYOU1", "HSP90B1", "CDC37", "SRC"],
  "title": "…",
  "bars": [{"bar": 1, "node": "HYOU1", "text": "…"}, "… 16 items, bars 1-4 node 1, 5-8 node 2, …"],
  "model": "claude-opus-5-5", "backend": "cli", "created": "2026-09-27T…Z"
}
```

Also in that folder: `walk_facts.json` (the data), `lyrics_prompt.md` (the prompt), `lyrics_raw.txt`, `lyrics.md`.
The lyrics are a creative summary of the data, not a scientific claim; a new run gives new lyrics.
