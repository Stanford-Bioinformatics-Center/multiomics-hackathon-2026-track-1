"""Tests for exvideo on a tiny made-up network (no pipeline outputs needed): python3 -m unittest discover network/video/tests"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from exvideo import align, audio, beats, dancer, lyrics, walk  # noqa: E402
from exvideo.errors import InputError, ModelError, WalkError  # noqa: E402
from exvideo.network import Edge, Network  # noqa: E402


def toy() -> Network:
    """A-B-C-D in a line, plus a strong B-E branch and an isolated-looking dead end E."""
    nodes = {n: {"node": n} for n in "ABCDE"}
    es = [("A", "B", 0.5, 0.1), ("B", "C", 0.2, 0.4), ("C", "D", 0.3, 0.3), ("B", "E", 0.9, 0.9)]
    return Network(nodes, {frozenset((a, b)): Edge(a, b, "protein - protein", x, y) for a, b, x, y in es})


def answer(walk_: list[str]) -> str:
    """A well-formed model answer for a walk."""
    return json.dumps({"title": "T", "suno_style": "g-funk", "sections": [{"node": n, "persona": "p", "bars": [f"{n} {i}" for i in range(4)]} for n in walk_]})


class WalkTests(unittest.TestCase):
    def test_resolve(self) -> None:
        net = toy()
        self.assertEqual(walk.resolve_node(" a ", net), "A")
        with self.assertRaises(WalkError):
            walk.resolve_node("Z", net)

    def test_check_walk(self) -> None:
        net = toy()
        walk.check_walk(["A", "B", "C", "D"], net)
        with self.assertRaisesRegex(WalkError, "not an edge"):
            walk.check_walk(["A", "C"], net)                       # no physical link: the hard gate
        with self.assertRaisesRegex(WalkError, "twice"):
            walk.check_walk(["A", "B", "A"], net)

    def test_random_walk(self) -> None:
        net = toy()
        w1, s1 = walk.random_walk("A", 3, net, seed=7)
        w2, _ = walk.random_walk("A", 3, net, seed=7)
        self.assertEqual(w1, w2)                                    # same seed, same walk
        self.assertEqual(s1, 7)
        self.assertEqual(w1, ("A", "B", "C", "D"))                  # the only 3-step path from A without revisits
        walk.check_walk(w1, net)
        with self.assertRaises(WalkError):
            walk.random_walk("E", 5, net, seed=1)                   # too few nodes for 5 steps


class LyricsTests(unittest.TestCase):
    def test_parse_ok(self) -> None:
        w = ["A", "B", "C", "D"]
        lyr = lyrics.parse_lyrics("here you go " + answer(w), w)
        self.assertEqual(len(lyr["bars"]), 16)
        self.assertEqual([b["bar"] for b in lyr["bars"]], list(range(1, 17)))
        self.assertEqual(lyr["bars"][4]["node"], "B")

    def test_parse_rejects(self) -> None:
        w = ["A", "B", "C", "D"]
        with self.assertRaisesRegex(ModelError, "should be"):
            lyrics.parse_lyrics(answer(["B", "A", "C", "D"]), w)    # wrong order
        bad = json.loads(answer(w)); bad["sections"][2]["bars"] = ["x", "y", "z"]
        with self.assertRaisesRegex(ModelError, "exactly 4"):
            lyrics.parse_lyrics(json.dumps(bad), w)
        with self.assertRaises(ModelError):
            lyrics.parse_lyrics("no json here", w)

    def test_prompt_is_verbatim(self) -> None:
        p = lyrics.build_prompt(["HYOU1", "HSP90B1", "CDC37", "SRC"], {"x": 1})
        self.assertTrue(p.startswith("give me 16 bars of 2pac rap lyrics which summarize the most interesting story from these 4 nodes: HYOU1 -> HSP90B1 -> CDC37 -> SRC."))
        self.assertIn("I. HYOU1: the first responder", p)           # the style example travels with it (neutral roles)

    def test_suno_text(self) -> None:
        w = ["A", "B", "C", "D"]
        s = lyrics.suno_text(lyrics.parse_lyrics(answer(w), w), w)
        self.assertIn("[Verse 1: A]", s)                          # lean tags (no persona text for Suno to sing)
        self.assertEqual(s.count("[Verse"), 4)


class AudioTests(unittest.TestCase):
    def test_timeline(self) -> None:
        w = ["A", "B", "C", "D"]
        bars = lyrics.parse_lyrics(answer(w), w)["bars"]
        tl = audio.timeline(bars, w, 60.0)
        self.assertAlmostEqual(tl["bars"][0]["start"], 6.0)
        self.assertAlmostEqual(tl["bars"][-1]["end"], 54.0)
        self.assertEqual([s["node"] for s in tl["segments"]], w)
        self.assertTrue(all(a["end"] <= b["start"] + 1e-9 for a, b in zip(tl["bars"], tl["bars"][1:])))
        with self.assertRaises(InputError):
            audio.timeline(bars, w, 10.0, intro=4, outro=4)

    def test_timeline_sung_headers(self) -> None:
        w = ["A", "B", "C", "D"]
        bars = lyrics.parse_lyrics(answer(w), w)["bars"]
        tl = audio.timeline(bars, w, 60.0, headers={n: f"{n} header" for n in w})
        self.assertEqual(len(tl["bars"]), 20)                        # 16 bars + 4 sung headers
        self.assertTrue(tl["bars"][0]["header"] and tl["bars"][5]["header"])
        self.assertEqual(tl["segments"][1]["start"], tl["bars"][5]["start"])   # a node's segment starts with its header

    def test_clean_path(self) -> None:
        self.assertEqual(audio.clean_path("'/tmp/my song.mp3' "), Path("/tmp/my song.mp3"))
        self.assertEqual(audio.clean_path("/tmp/my\\ song.mp3"), Path("/tmp/my song.mp3"))


class AlignTests(unittest.TestCase):
    def test_normalise(self) -> None:
        self.assertEqual(align.normalise("Four hours, I'm the ONE"), ["4", "hours", "im", "the", "1"])

    def test_sync_lines(self) -> None:
        lines = [{"node": "A", "text": "I. A: the lookout", "header": True}, {"node": "A", "text": "four hours after the whistle"},
                 {"node": "A", "text": "endurance or the iron"}, {"node": "B", "text": "six sugars on my jacket"}]
        heard = "4 hours after the whistle endurance of the eye and 6 sugars on my jacket".split()
        words = [{"word": w, "start": 10.0 + i, "end": 10.8 + i} for i, w in enumerate(heard)]   # mishearings on purpose
        timed, rep = align.sync_lines(lines, words, 40.0)
        self.assertEqual(rep["dropped_headers"], ["I. A: the lookout"])              # not sung: left out
        self.assertEqual([t["start"] for t in timed], [10.0, 15.0, 20.0])            # each line starts at its first word ("6" is the 11th heard word)
        self.assertEqual(timed[0]["end"], timed[1]["start"])
        tl = audio.timeline_from_lines(timed, ["A", "B"], 40.0)
        self.assertEqual(tl["intro"], 10.0)
        self.assertEqual(tl["segments"][0]["end"], tl["segments"][1]["start"])


class DancerTests(unittest.TestCase):
    """A synthetic bobbing block: two colour clips with a solid flash frame between; hits must be the bottoms of the bob."""

    @staticmethod
    def _gif():
        import numpy as np
        frames = []
        for i in range(45):
            f = np.zeros((60, 30, 4), np.uint8); y = int(round(20 + 8 * np.sin(2 * np.pi * (i - 2.5) / 10)))   # lowest at 5, 15, 25, ...
            f[y:y + 20, 10:20] = (120, 120, 110, 255) if i < 30 else (130, 70, 65, 255)
            frames.append(f)
        flash = np.full((60, 30, 4), (0, 200, 0, 255), np.uint8)
        return np.stack(frames[:30] + [flash] + frames[30:])

    def test_clips_drop_flash_and_split_on_colour(self):
        g = dancer.clips(self._gif())
        self.assertEqual([len(c) for c in g], [30, 15])
        self.assertNotIn(30, g[0] + g[1])                  # the flash frame is gone

    def test_hits_and_seamless_even_loop(self):
        a = self._gif()[:30, ..., 3] > 128
        h = dancer.hits(a)
        self.assertEqual(h, [5, 15, 25])
        self.assertEqual(dancer.best_loop(a, h), (0, 2))   # two steps, end pose identical to the start pose

    def test_despill_removes_green_fringe_only(self):
        import numpy as np
        px = np.array([[[90, 160, 80, 255], [120, 120, 110, 255], [0, 255, 0, 0]]], np.uint8)
        out = dancer.despill(px)
        self.assertEqual(out[0, 0].tolist(), [90, 90, 80, 255])
        self.assertEqual(out[0, 1].tolist(), [120, 120, 110, 255])
        self.assertEqual(out[0, 2].tolist(), [0, 255, 0, 0])   # fully transparent pixels untouched

class BeatTests(unittest.TestCase):
    """Synthetic drum grooves (kick every beat, snare on 2 and 4, hi-hat on the off-beats): the grooves that fooled the
    old autocorrelation tempo (half time, 1.5-beat lags) and the old whole-spectrum flux (locked onto the hi-hats)."""

    @staticmethod
    def _groove(bpm, dur=30.0, seed=1):
        import numpy as np
        sr, rng = beats.SR, np.random.default_rng(seed)
        y = np.zeros(int(dur * sr)); tk = np.arange(int(0.25 * sr)) / sr
        kick = np.sin(2 * np.pi * (50 + 80 * np.exp(-tk * 30)) * tk) * np.exp(-tk * 12)
        hat = rng.standard_normal(int(0.04 * sr)) * np.exp(-np.arange(int(0.04 * sr)) / sr * 120) * 0.3
        snare = rng.standard_normal(int(0.15 * sr)) * np.exp(-np.arange(int(0.15 * sr)) / sr * 25) * 0.5
        true = np.arange(0.5, dur - 1, 60 / bpm)
        for k, b in enumerate(true):
            for sound, t in ((kick, b), (snare, b if k % 2 else None), (hat, b + 30 / bpm)):
                if t is not None:
                    i = int(t * sr); y[i:i + len(sound)] += sound[:max(0, len(y) - i)]
        return (y / np.abs(y).max()).astype(np.float32), true

    def _detect(self, y):
        import numpy as np
        o = beats.onset_strength(y); p = beats.tempo_period(o)
        return 60 * beats.SR / beats.HOP / p, (beats.track(o, p) * beats.HOP + beats.NFFT / 2) / beats.SR + beats.ONSET_LAG

    def test_tempo_and_phase_across_tempi(self):
        import numpy as np
        for bpm in (92, 128, 174):
            y, true = self._groove(bpm)
            est, got = self._detect(y)
            self.assertLess(abs(est - bpm) / bpm, 0.02, f"{bpm} BPM detected as {est:.1f}")
            err = [np.min(np.abs(got - b)) for b in true[2:-2]]
            self.assertLess(max(err), 0.03, f"{bpm} BPM: a beat is {1000 * max(err):.0f} ms off")      # on the beat, not the off-beat

    def test_frame_at_puts_hits_on_beats(self):
        grid, hits, n = [0.0, 0.6, 1.3, 1.9, 2.6], [0, 17], 35
        self.assertEqual([dancer.frame_at(t, grid, hits, n, 1.0, 30) for t in grid[:4]], [0, 17, 0, 17])
        self.assertEqual(dancer.frame_at(0.3, grid, hits, n, 1.0, 30), 8)                 # half way between two hits
        self.assertEqual(dancer.frame_at(0.6 - 0.01, grid, hits, n, 1.0, 30), 17)         # nearest video frame snaps to the hit
        self.assertEqual([dancer.frame_at(t, grid, hits, n, 0.5, 30) for t in grid[:3]], [0, 8, 17])   # a step every other beat

class CostumeTests(unittest.TestCase):
    """A synthetic grey head with two pink ears and a pink nose: the band goes under the ears, the knot above the nose."""

    def test_band_under_ears_knot_above_nose(self):
        import numpy as np
        from exvideo import costume
        f = np.zeros((60, 40, 4), np.uint8)
        f[10:50, 8:32] = (120, 120, 110, 255)                          # head + body
        f[10:14, 10:14] = (200, 120, 120, 255); f[10:14, 26:30] = (200, 120, 120, 255)   # ears
        f[24:27, 23:26] = (200, 120, 120, 255)                          # nose, right of centre
        g = costume.head_geometry(f)
        self.assertTrue(13 < g["ly"] < 17 and 13 < g["ry"] < 17, g)     # just under the ears
        self.assertAlmostEqual(g["nose"], 24.0, delta=0.6)
        out = costume.bandana(np.stack([f] * 4))
        self.assertEqual(out.shape, (4, 60 * costume.UP, 40 * costume.UP, 4))
        cloth = (out[0, ..., 0] > 70) & (out[0, ..., 1] < 70) & (out[0, ..., 3] > 200)      # red cloth, not the pink ears
        top = (min(g["ly"], g["ry"]) - g["thick"] / 2 - 1.0) * costume.UP
        ys, xs = np.nonzero(cloth); ends = xs[ys < top]                  # only the knot's ends rise above the band
        self.assertTrue(np.any(cloth[int(g["ly"] * costume.UP)]))       # the band is there
        self.assertGreater(len(ends), 0)
        self.assertAlmostEqual(ends.mean() / costume.UP, 24.0, delta=3.0)   # the knot is above the nose (2Pac: in front)

class SunoBanTests(unittest.TestCase):
    """Suno rejects artist / producer names and producer tags ("producer tag phosphate"): the list is in the prompt,
    and every answer is checked; only the offending bars are rewritten."""

    def _lyr(self, texts):
        return {"title": "T", "suno_style": "g-funk", "personas": {"A": "the boss"}, "bars": [{"bar": i + 1, "node": "A", "text": t} for i, t in enumerate(texts)]}

    def test_find_banned_whole_words_and_plurals(self):
        hits = lyrics.find_banned(self._lyr(["Ten phosphates on file", "my phosphatase and pacman", "Snoop on the line"]), ["phosphate", "snoop", "pac"])
        self.assertEqual(hits, [("bar 1", "phosphates"), ("bar 3", "Snoop")])

    def test_prompt_carries_the_rules_and_list(self):
        p = lyrics.build_prompt(["A", "B", "C", "D"], {})
        self.assertIn("SUNO RULES", p); self.assertIn("phosphate", p.split("SUNO RULES")[1])
        self.assertTrue(p.startswith(lyrics.PROMPT_TEMPLATE.format(walk="A -> B -> C -> D")))   # the team's words stay verbatim, first

    def test_scrub_rewrites_only_offending_bars(self):
        walk = ["A"]; old = lyrics.BACKENDS.get("cli")
        answer = json.dumps({"title": "NEW", "suno_style": "x", "sections": [{"node": "A", "persona": "p", "bars": [
            "changed line one that nobody should ever keep today", "Ten P-tags on the file and every one is on the clock,",
            "changed line three that nobody should ever keep today", "changed line four that nobody should ever keep today"]}]})
        saved = lyrics.BARS_PER_NODE
        try:
            lyrics.BACKENDS["cli"] = lambda prompt, model: answer
            lyr = self._lyr(["Four hours after the whistle I'm the first one on the block,", "Ten phosphate spots on file and every one is on the clock,",
                             "Endurance or the iron, either way my numbers climb,", "Only three lines on my phone but I'm the spark every time,"])
            out = lyrics.scrub(lyr, walk, "cli", "m")
        finally:
            lyrics.BACKENDS["cli"] = old
        self.assertEqual([b["text"] for b in out["bars"]], ["Four hours after the whistle I'm the first one on the block,", "Ten P-tags on the file and every one is on the clock,",
                                                      "Endurance or the iron, either way my numbers climb,", "Only three lines on my phone but I'm the spark every time,"])
        self.assertEqual(out["title"], "T")

class LengthTests(unittest.TestCase):
    """Songs must stay under 1:30: bars of at most MAX_WORDS words, short personas; the prompt says so."""

    def test_find_long(self):
        lyr = {"personas": {"A": "the runner who moves the lactic out"}, "bars": [{"bar": 1, "node": "A", "text": "short line here"},
               {"bar": 2, "node": "A", "text": " ".join(["word"] * (lyrics.MAX_WORDS + 1))}]}
        self.assertEqual(lyrics.find_long(lyr), [("bar 2", f"{lyrics.MAX_WORDS + 1} words"), ("persona A", "7 words")])

    def test_prompt_has_length_rules(self):
        p = lyrics.build_prompt(["A", "B", "C", "D"], {})
        self.assertIn("LENGTH RULES", p); self.assertIn(f"at most {lyrics.MAX_WORDS} words", p)

class RespectAndRhymeTests(unittest.TestCase):
    """No gang / street-crime words (personas, bars, style); bars 10-16 syllables and AABB on the last word's vowel."""

    def test_respect_words_caught(self):
        lyr = {"title": "T", "suno_style": "90s gangsta rap", "personas": {"A": "the kingpin"}, "bars": [{"bar": 1, "node": "A", "text": "I run the block like a hub"}]}
        where = [w for w, _ in lyrics.find_banned(lyr, lyrics.banned_terms("respect"))]
        self.assertEqual(where, ["style", "persona A"])

    def test_rhyme_and_syllables(self):
        good = ["Four hours after the whistle I'm the first one on the block,", "Endurance or the iron, either way my numbers pop,",
                "Six sugars on my jacket, N-linked, I stay dressed,", "Diabetes checked my papers and I passed the test,"]
        bad = ["Thirty-four sugars on my coat and twenty-three are N-linked", "Years stack me higher but more of me keeps diabetes low",
               "Weights hit and by four hours I dipped inside the blood", "One blue line"]
        mk = lambda t: {"bars": [{"bar": i + 1, "node": "A", "text": x} for i, x in enumerate(t)]}
        self.assertEqual(lyrics.find_rhythm(mk(good), ["A"]), [])
        why = lyrics.find_rhythm(mk(bad), ["A"])
        self.assertIn(("bar 4", "3 syllables"), why)
        self.assertTrue(any("does not rhyme" in t for _, t in why))

if __name__ == "__main__":
    unittest.main()
