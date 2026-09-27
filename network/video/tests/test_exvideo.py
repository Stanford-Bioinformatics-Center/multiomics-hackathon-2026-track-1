"""Tests for exvideo on a tiny made-up network (no pipeline outputs needed): python3 -m unittest discover network/video/tests"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from exvideo import align, audio, dancer, lyrics, walk  # noqa: E402
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
        self.assertIn("I. HYOU1: the lookout", p)                    # the style example travels with it

    def test_suno_text(self) -> None:
        w = ["A", "B", "C", "D"]
        s = lyrics.suno_text(lyrics.parse_lyrics(answer(w), w), w)
        self.assertIn("[Verse 1: I. A, p]", s)
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

if __name__ == "__main__":
    unittest.main()
