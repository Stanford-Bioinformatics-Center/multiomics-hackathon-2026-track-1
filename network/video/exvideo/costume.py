"""The rat's costume: a red Team 2-PAC bandana tied round its head, drawn onto every frame of the dance loop.

How it works (plain language):
  1. every frame is enlarged 4x first (Lanczos, alpha premultiplied so the cut-out edge does not darken), because the
     rat's head is only ~20 pixels wide in the GIF; the bandana is then drawn at that size, 2x supersampled (smooth);
  2. the head is found frame by frame from the rat itself: the top of the silhouette, the two pink EARS in the top rows
     and the pink NOSE below them;
  3. the band runs across the forehead just under the ears, tilted like the line between the ears, and is clipped to the
     head's outline (so it wraps round the head rather than floating); the KNOT is tied 2Pac-style at the FRONT, on
     the forehead right above the nose, with its two ends sticking up and out (they swing a little against the rat's
     movement);
  4. all of that geometry is smoothed over neighbouring frames (the loop wraps) so the bandana does not jitter;
  5. the cloth is cardinal red with darker folds at the edges, a highlight and a small cream paisley-dot print that is
     fixed to the band, so it moves with it.
The colours follow the team's Gemini badge art (Downloads, 2026-09-27): a cardinal paisley bandana; the tie is 2Pac's.
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

from .errors import InputError

UP = 4                          # enlargement of the GIF before drawing
SS = 2                          # supersampling of the bandana layer
CLOTH = (140, 21, 21)           # cardinal red (#8C1515)
FOLD = (84, 10, 12)             # the shaded folds / edges
SHINE = (196, 58, 52)           # the highlight along the band
PRINT = (244, 232, 208)         # the cream paisley print


def _pil():
    try:
        from PIL import Image, ImageDraw, ImageFilter
    except ImportError:
        raise InputError("the bandana needs Pillow: network/video/.venv/bin/python -m pip install pillow (or bash network/video/setup.sh)")
    return Image, ImageDraw, ImageFilter


def _clusters(xs: np.ndarray, ys: np.ndarray, gap: int = 3) -> List[Tuple[float, float, int, int]]:
    """Pink pixels grouped by x (a gap of `gap` columns splits groups): (centre x, lowest y, size, top y) per group."""
    if len(xs) == 0:
        return []
    order = np.argsort(xs); xs, ys = xs[order], ys[order]
    cut = np.where(np.diff(xs) > gap)[0] + 1
    return [(float(gx.mean()), float(gy.max()), len(gx), int(gy.min())) for gx, gy in zip(np.split(xs, cut), np.split(ys, cut))]


def head_geometry(rgba: np.ndarray) -> Dict[str, float]:
    """Where the band goes on one frame (GIF pixel units): left/right end points of its centre line, its thickness, and
    the nose's x (the front of the face, where the knot is tied)."""
    op = rgba[..., 3] > 128
    rows = np.where(op.any(axis=1))[0]
    if len(rows) == 0:
        raise InputError("empty dancer frame")
    y0 = int(rows[0])
    r, g, b = (rgba[..., k].astype(int) for k in range(3))
    pink = op & (r - g > 25) & (r - b > 25)
    ys, xs = np.nonzero(pink[y0:y0 + 10]); ears = sorted(_clusters(xs, ys + y0), key=lambda c: -c[2])[:2]
    ys, xs = np.nonzero(pink[y0 + 10:y0 + 24]); nose = sorted(_clusters(xs, ys + y0 + 10), key=lambda c: -c[2])[:1]
    head = op[y0:y0 + 22]; cols = np.nonzero(head.any(axis=0))[0]
    hx = float(cols.mean()) if len(cols) else rgba.shape[1] / 2
    if len(ears) == 2:
        (xa, ya), (xb, yb) = sorted(((e[0], e[1]) for e in ears))
    elif len(ears) == 1:
        xa = xb = ears[0][0]; ya = yb = ears[0][1]
    else:
        xa = xb = hx; ya = yb = y0 + 5.0
    ya += 1.8; yb += 1.8                                    # just under the ears
    slope = (yb - ya) / (xb - xa) if xb - xa > 2 else 0.0
    slope = float(np.clip(slope, -0.6, 0.6))
    ymid = (ya + yb) / 2; xmid = (xa + xb) / 2
    row = op[int(round(ymid))]; xs_row = np.nonzero(row)[0]
    if len(xs_row):                                          # the band spans the head's width at that height, + 1 px each side
        lx, rx = float(xs_row[xs_row < xmid + 12].min()) - 1.0, float(xs_row[xs_row > xmid - 12].max()) + 1.0
    else:
        lx, rx = xmid - 8.0, xmid + 8.0
    return {"lx": lx, "ly": ymid + slope * (lx - xmid), "rx": rx, "ry": ymid + slope * (rx - xmid),
            "thick": max(2.6, 0.2 * (rx - lx)), "nose": nose[0][0] if nose else hx, "cx": float(np.nonzero(op)[1].mean())}


def _smooth(geo: List[Dict[str, float]], k: int = 2) -> List[Dict[str, float]]:
    """Circular moving average (the loop wraps) of every number, over 2k+1 frames."""
    n = len(geo)
    out = [{key: float(np.mean([geo[(i + d) % n][key] for d in range(-k, k + 1)])) for key in geo[i]} for i in range(n)]
    for i in range(n):                                       # how fast the body moves sideways (for the tails' swing)
        out[i]["vx"] = (geo[(i + 1) % n]["cx"] - geo[i - 1]["cx"]) / 2
    return out


def _upscale(frame: np.ndarray):
    Image, _, _ = _pil()
    h, w = frame.shape[:2]
    return Image.fromarray(frame, "RGBA").convert("RGBa").resize((w * UP, h * UP), Image.LANCZOS).convert("RGBA")


def _draw(img, g: Dict[str, float], i: int):
    """The bandana on one enlarged frame."""
    Image, ImageDraw, ImageFilter = _pil()
    W, H = img.size; S = UP * SS
    layer = Image.new("RGBA", (W * SS, H * SS), (0, 0, 0, 0)); d = ImageDraw.Draw(layer)
    P = lambda x, y: (x * S, y * S)                                           # GIF pixel units -> layer pixels
    lx, ly, rx, ry, th = g["lx"], g["ly"], g["rx"], g["ry"], g["thick"]
    ux, uy = rx - lx, ry - ly; L = max(np.hypot(ux, uy), 1e-6); ux, uy = ux / L, uy / L; nx, ny = uy, -ux     # along, up
    bulge = 0.9                                                               # the forehead is round: the band arches up
    def edge(off, m=24):
        return [P(lx + ux * L * t + nx * (off + bulge * 4 * t * (1 - t)), ly + uy * L * t + ny * (off + bulge * 4 * t * (1 - t))) for t in np.linspace(0, 1, m)]
    top, bot = edge(th / 2), edge(-th / 2)
    band = Image.new("L", layer.size, 0); ImageDraw.Draw(band).polygon(top + bot[::-1], fill=255)
    d.polygon(top + bot[::-1], fill=CLOTH + (255,))
    d.line(edge(th * 0.12), fill=SHINE + (255,), width=max(1, int(S * th * 0.22)))                    # highlight
    d.line(top, fill=FOLD + (255,), width=max(1, int(S * 0.55))); d.line(bot, fill=FOLD + (255,), width=max(1, int(S * 0.7)))
    for t in np.arange(0.08, 0.95, 0.14):                                     # the paisley print: a dot and a comma
        off = (0.18 if int(t * 50) % 2 else -0.18) * th
        x, y = lx + ux * L * t + nx * (off + bulge * 4 * t * (1 - t)), ly + uy * L * t + ny * (off + bulge * 4 * t * (1 - t))
        rr = 0.33 * S
        d.ellipse([x * S - rr, y * S - rr, x * S + rr, y * S + rr], fill=PRINT + (235,))
        d.arc([x * S + 0.2 * S, y * S - 0.9 * S, x * S + 1.5 * S, y * S + 0.4 * S], 200, 340, fill=PRINT + (200,), width=max(1, int(0.28 * S)))
    # clip the band to the head (dilated 1 px) so it wraps round it
    alpha = img.split()[3].resize(layer.size, Image.BILINEAR).point(lambda a: 255 if a > 90 else 0).filter(ImageFilter.MaxFilter(2 * S + 1))
    from PIL import ImageChops
    layer.putalpha(ImageChops.multiply(layer.split()[3], ImageChops.multiply(band.filter(ImageFilter.MaxFilter(3)), alpha)))
    # the knot, 2Pac-style: in FRONT, on the band right above the nose, its two ends sticking up and out
    t = float(np.clip((g["nose"] - lx) / max(rx - lx, 1e-6), 0.3, 0.7)); lift = bulge * 4 * t * (1 - t)
    kx, ky = lx + ux * L * t + nx * lift, ly + uy * L * t + ny * lift
    knot = Image.new("RGBA", layer.size, (0, 0, 0, 0)); k = ImageDraw.Draw(knot)
    swing = float(np.clip(-g["vx"] * 14, -25, 25)); tilt = np.degrees(np.arctan2(uy, ux))
    for j, (ang, ln, wd) in enumerate(((-145, 6.4, 2.9), (-35, 6.8, 2.9))):   # up-left and up-right of the band
        a = np.radians(ang + tilt + swing + 5 * np.sin(2 * np.pi * i / 12 + 2 * j))
        dx, dy = np.cos(a), np.sin(a); px, py = -dy, dx                         # along / across the end
        tip = (kx + dx * ln, ky + dy * ln)
        end = [P(kx + px * wd * 0.35, ky + py * wd * 0.35), P(tip[0] + px * wd * 0.55, tip[1] + py * wd * 0.55),
               P(tip[0] + dx * 0.6, tip[1] + dy * 0.6), P(tip[0] - px * wd * 0.45, tip[1] - py * wd * 0.45),
               P(kx - px * wd * 0.35, ky - py * wd * 0.35)]
        k.polygon(end, fill=CLOTH + (255,), outline=FOLD + (255,))
        k.line([P(kx, ky), P(tip[0] - px * wd * 0.05, tip[1] - py * wd * 0.05)], fill=FOLD + (200,), width=max(1, int(0.35 * S)))   # the crease
    r0 = th * 0.62 * S
    k.ellipse([kx * S - r0, ky * S - r0 * 0.85, kx * S + r0, ky * S + r0 * 0.85], fill=CLOTH + (255,), outline=FOLD + (255,), width=max(1, int(0.5 * S)))
    k.ellipse([kx * S - r0 * 0.35, ky * S - r0 * 0.5, kx * S + r0 * 0.15, ky * S - r0 * 0.05], fill=SHINE + (255,))
    layer = Image.alpha_composite(layer, knot)
    return Image.alpha_composite(img, layer.resize((W, H), Image.LANCZOS))


def bandana(loop: np.ndarray) -> np.ndarray:
    """The dance loop (frames, h, w, 4) -> the same loop enlarged UP x, each frame wearing the bandana."""
    geo = _smooth([head_geometry(f) for f in loop])
    return np.stack([np.array(_draw(_upscale(f), g, i)) for i, (f, g) in enumerate(zip(loop, geo))])
