"""Tests for exvideo on a tiny made-up network (no pipeline outputs needed): python3 -m unittest discover network/video/tests"""
from __future__ import annotations

import json
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from exvideo import audio, lyrics, walk  # noqa: E402
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


if __name__ == "__main__":
    unittest.main()
