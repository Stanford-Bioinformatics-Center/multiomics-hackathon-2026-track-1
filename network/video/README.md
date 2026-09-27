# Music video (Team 2-PAC)

**One command, about five minutes including Suno:** pick a node → random walk through the network → 16 bars of lyrics
written from the figure 17 data → you make the song in Suno → a fly-through video of the walk with the lyrics on screen.

```bash
python3 network/video/make_music_video.py            # asks for a node, writes the lyrics, waits for the song, renders
python3 network/video/make_music_video.py --start HYOU1 --seed 7                 # reproducible random walk
python3 network/video/make_music_video.py --walk HYOU1,HSP90B1,CDC37,SRC         # a fixed walk
python3 network/video/make_music_video.py --start SRC --audio ~/Downloads/song.mp3   # song already made
```

**In two sittings** (the song is made later, or by someone else):

```bash
python3 network/video/01_lyrics_from_walk.py          # asks for a node -> the team's random walk -> 16 bars -> clipboard; stops
python3 network/video/02_video_from_song.py <the walk folder it printed> ~/Downloads/song.mp3     # the song back -> the video
```
The start node must be a node of figure 17 (the 17a joint network) that a full 4-node walk can start from (315 of the
353; the prompt suggests close names and refuses dead ends; Enter = a random start). The walk is the team's
`random_walk/random_walks.R` with a **coin flip at every step** (`--arm coin`, the default): heads endurance, tails
resistance weights choose that step, so the likely paths stay likely but the arms mix. It only takes steps from which
the walk can still be finished, so every walk has exactly 4 different nodes (no dead ends, no repeats, no retries).
`--arm EE` / `--arm RE` use one arm throughout. Part 1 takes every `make_music_video.py` option (`--start HYOU1 --seed 7`, `--walk ...`); it saves everything
part 2 needs in `$HACK_OUT/video/<walk>/` and prints the exact part-2 command. Leave the song out of part 2 and it asks
for it (drag it in) or waits for the next download.

**Suno rejected a word?** ("Your lyrics contain producer tag phosphate - we don't reference specific artists")
`python3 network/video/fix_lyrics.py <walk folder> phosphate` bans the word from now on (`exvideo/suno_banned.txt`)
and rewrites only the bars that use it; everything else stays word for word, and the fixed lyrics are printed and put
on the clipboard. Every new lyric run gets the banned list in its prompt (after the team's words, which stay
verbatim) and is checked against it before it is shown; offending bars are rewritten automatically.

**Just the rat, on any song:** `rat_dance.py` makes a video of the dancing rat stepping exactly on the beats of whatever
song you give it (about 7 s for a 90 s song; no Node, no network data):

```bash
python3 network/video/rat_dance.py ~/Downloads/song.mp3               # -> ~/Downloads/song_rat.mp4 (720x1280, white)
python3 network/video/rat_dance.py song.mp3 --size 1920x1080 --bg black --out ~/Desktop/rat.mp4
python3 network/video/rat_dance.py song.mp3 --transparent             # ProRes 4444 .mov with alpha, to drop into an edit
python3 network/video/rat_dance.py song.mp3 --bpm 170                 # if it dances at half / double time
python3 network/video/rat_dance.py song.mp3 --nudge -0.03             # steps 30 ms earlier
python3 network/video/rat_dance.py song.mp3 --no-bandana              # the plain rat
```
The rat wears the red Team 2-PAC bandana, tied 2Pac-style with the knot in front (`exvideo/costume.py`: the GIF is
enlarged 4x; the band is fitted once, under the ears of the most typical frame, and then carried by the head's own
motion, found by matching the whole head frame to frame, so it stays put; clipped to the head; knot above the nose with
both ends up and out; colours from the team's badge art). Same rat, same flag,
in the music video.
Checked on synthetic drum grooves at 70, 92, 128 and 174 BPM and a 100→120 BPM ramp: tempo found in all five (70 at
double time), every beat within 14 ms, and the rat on a step-hit pose at every beat it steps on.

## What happens

| # | Stage | Time |
|---|---|---|
| 1 | You type a node / feature (gene symbol or metabolite; close matches are suggested) | — |
| 2 | **Random walk**, 3 steps, with the **team's walker** (`random_walk/random_walks.R`, Subarna Bhattacharya; see its README): each step goes to a not-yet-visited neighbour with probability set by the chosen arm's edge weights (`--arm EE` endurance, the default, or `RE`); dead ends are redrawn; seeded and saved with the step probabilities. Every step is re-checked against the physical edges. (`--walker builtin`: exvideo's own walker, probability ∝ max(\|w_EE\|, \|w_RE\|)) | ~1 s |
| 3 | **Lyrics**: the figure 17 facts of the 4 nodes and 3 edges (hub rank, strength per arm, module pathway, strongest exercise response, T2D and ageing directions from Öhman 2021 and UK Biobank, phosphosites, sugars, top partners; per edge the weights, their difference and arm specificity) + **the team's prompt, verbatim** + the team's favourite lyrics as the **style example** → Claude (Claude Code `claude -p`, or `--backend api`) → checked: 4 bars per node in walk order, one retry if malformed | ~25 s |
| 4 | **Suno**: the Suno-formatted lyrics are copied to the clipboard, the style prompt and title printed, suno.com/create opened. In Suno: *Create → Custom*, paste, set the style, create, download | ~1-2 min (you) |
| 5 | **The song**: drag the file into the terminal, or press Enter and it picks up the new download in `~/Downloads` | — |
| 5b | **Sync**: Whisper hears the song; every lyric line is placed where its first word is sung (first run ~45 s per 80 s of song, then cached) | ~45 s |
| 5c | **Beat**: tempo and every beat of the song (spectral-flux onsets in 40 log-spaced bands, pulse-train tempo over 72-176 BPM, dynamic-programming beat tracking; numpy); the dancer GIF cleaned (flash frames dropped, one single-colour clip kept, green fringe removed) and its step hits found (the bottom of each bob) with a seamless even-step loop | ~1 s |
| 6 | **The video** (`render/render_walk.js`, Puppeteer + ffmpeg): the figure 17 page full screen with arm-specific edges on (red = endurance only, blue = resistance only); title card; the camera flies node to node, each walked edge turns gold; each node's persona and fact card; the current bar large and the next bar faded; the whole walk at the end; a **dancing rat** at the side played hit to hit: each step's hit frame is shown exactly on a beat and the frames in between are spread over the beat, so it follows the song and never drifts (grey clip by default, `--dancer-clip 1` red / `2` teal; `--dancer-steps-per-beat 0.5|1|2`; credit only on the closing card; `--dancer other.gif`, `--no-dancer`); the song underneath → `music_video.mp4` | ~1 s per second of song |

**Lyrics are synced to the actual vocals** (`--sync whisper`, the default): Whisper (faster-whisper, run locally)
transcribes the song with a time for every word; the heard words are then aligned to the known lyrics (a global word
alignment that tolerates mishearings such as "HU1" for HYOU1 or "947" for "nine-forty-seven"), and every line starts
when its first word is sung. Section headers that are not sung are left out automatically; the beat before the first
line shows the title card and the beat after the last line the full walk. The transcript is cached next to the video,
so re-renders skip it. `sync_report.json` shows how many lyric words were heard. `--sync even` spreads the lines evenly
instead (no Whisper). To film a song you already have: `--lyrics <lyrics.json> --audio <song>`.
Suno style prompts must not name artists; the one Claude writes does not.

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
`python3 -m unittest discover -s network/video/tests` (17 tests on a toy network and synthetic drum grooves, incl. lyric-to-vocal alignment, beat tracking and the dancer: walks, the hard gate, reproducible
random walks, answer checks, the verbatim prompt, Suno text, timing, dragged paths).

**The dancer GIF** (not in the repo): the rat-dance meme (original by @ratomilton, TikTok), transparent GIF from Tenor,
saved at `~/Desktop/output/hackathon-2026-track1/external/dancer/rat_dance_transparent.gif` (override with `--dancer`).

**One-time set-up:** `bash network/video/setup.sh` (Whisper environment in `network/video/.venv`, renderer packages;
`make_music_video.py` switches to `.venv` by itself).

**Needs:** the pipeline outputs incl. the figure 17 page (`bash network/run_all.sh`); Claude Code (`claude`) or
`ANTHROPIC_API_KEY`; Node.js (the renderer installs Puppeteer on first use); ffmpeg (`brew install ffmpeg`).
The lyrics are a creative summary of the data, not a scientific claim.
