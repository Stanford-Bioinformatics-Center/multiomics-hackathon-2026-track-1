# Music video (Team 2-PAC)

**One command, about five minutes including Suno:** pick a node → random walk through the network → 16 bars of lyrics
written from the figure 17 data → you make the song in Suno → a fly-through video of the walk with the lyrics on screen.

```bash
python3 network/video/make_music_video.py            # asks for a node, writes the lyrics, waits for the song, renders
python3 network/video/make_music_video.py --start HYOU1 --seed 7                 # reproducible random walk
python3 network/video/make_music_video.py --walk HYOU1,HSP90B1,CDC37,SRC         # a fixed walk
python3 network/video/make_music_video.py --start SRC --audio ~/Downloads/song.mp3   # song already made
python3 network/video/01_lyrics_from_walk.py --start CDC37                       # lyrics only (no Suno, no video)
```

## What happens

| # | Stage | Time |
|---|---|---|
| 1 | You type a node / feature (gene symbol or metabolite; close matches are suggested) | — |
| 2 | **Random walk**, 3 steps along the network's physical edges (STRING / Rhea): each step goes to a not-yet-visited neighbour with probability ∝ max(\|w_EE\|, \|w_RE\|), so strongly co-regulated links are likelier; seeded and saved | < 1 s |
| 3 | **Lyrics**: the figure 17 facts of the 4 nodes and 3 edges (hub rank, strength per arm, module pathway, strongest exercise response, T2D and ageing directions from Öhman 2021 and UK Biobank, phosphosites, sugars, top partners; per edge the weights, their difference and arm specificity) + **the team's prompt, verbatim** + the team's favourite lyrics as the **style example** → Claude (Claude Code `claude -p`, or `--backend api`) → checked: 4 bars per node in walk order, one retry if malformed | ~25 s |
| 4 | **Suno**: the Suno-formatted lyrics are copied to the clipboard, the style prompt and title printed, suno.com/create opened. In Suno: *Create → Custom*, paste, set the style, create, download | ~1-2 min (you) |
| 5 | **The song**: drag the file into the terminal, or press Enter and it picks up the new download in `~/Downloads` | — |
| 6 | **The video** (`render/render_walk.js`, Puppeteer + ffmpeg): the figure 17 page full screen with arm-specific edges on (red = endurance only, blue = resistance only); title card; the camera flies node to node, each walked edge turns gold; each node's persona and fact card; the current bar large and the next bar faded; the whole walk at the end; the song underneath → `music_video.mp4` | ~1 s per second of song |

If the song also sings the section headers ("I. HYOU1: the lookout"), add `--sung-headers`: each node's segment then
starts with a slot for its header (shown in gold). To film a song you already have: `--lyrics <lyrics.json> --audio <song>`
(skips the lyrics and Suno steps). Bars are spread evenly over the vocal part of the song (10% intro and outro, at most 8 s each; tune with
`--intro` / `--outro` and re-render). Suno style prompts must not name artists; the one Claude writes does not.

## Files (`$HACK_OUT/video/<walk>/`)

`walk.json` (walk + seed) · `walk_facts.json` (the data sent to the model) · `lyrics_prompt.md` (the exact prompt) ·
`lyrics_raw.txt` · **`lyrics.json`** (the contract below) · `lyrics.md` · `suno_lyrics.txt` · `suno_style.txt` ·
`render_spec.json` · **`music_video.mp4`**

```json
{"walk": ["HYOU1", "HSP90B1", "CDC37", "SRC"], "seed": 7, "title": "…", "suno_style": "…",
 "personas": {"HYOU1": "…", "…": "…"}, "audio": "…/song.mp3", "duration": 45.0, "intro": 4.5, "outro": 4.5,
 "bars": [{"bar": 1, "node": "HYOU1", "text": "…", "start": 4.5, "end": 6.75}, "… 16 bars"],
 "segments": [{"node": "HYOU1", "start": 4.5, "end": 13.5}, "… one per node"], "model": "claude-opus-5-5", "backend": "cli"}
```

## Code

`exvideo/` (typed Python package): `walk.py` (node lookup, walk checks — every step must be a physical edge —
random walk), `network.py` (the joint network, the figure 17 facts, the on-screen fact card), `lyrics.py` (the team's
prompt, the style example, the Claude call, answer checks, Suno formatting), `audio.py` (waiting for the download,
song length, bar timing), `render.py` (renderer set-up, frames → MP4), `errors.py` (classed errors).
`render/render_walk.js` draws the frames (frame-stepped, so none are dropped). Tests:
`python3 -m unittest discover -s network/video/tests` (9 tests on a toy network: walks, the hard gate, reproducible
random walks, answer checks, the verbatim prompt, Suno text, timing, dragged paths).

**Needs:** the pipeline outputs incl. the figure 17 page (`bash network/run_all.sh`); Claude Code (`claude`) or
`ANTHROPIC_API_KEY`; Node.js (the renderer installs Puppeteer on first use); ffmpeg (`brew install ffmpeg`).
The lyrics are a creative summary of the data, not a scientific claim.
