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
  4. the band is fitted ONCE, on the most typical frame, and then carried by the head's own rigid motion (rotation +
     shift found by matching the whole head against that frame), so it stays put on the head and does not wobble;
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
TILT_KEEP, TILT_MAX = 0.4, 8.0  # share of the head's fitted rotation the band follows, and its cap (degrees)


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
    t = g["knot_t"] if "knot_t" in g else float(np.clip((g["nose"] - lx) / max(rx - lx, 1e-6), 0.3, 0.7)); lift = bulge * 4 * t * (1 - t)
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


def _features(frame: np.ndarray) -> np.ndarray:
    """What the head tracker compares: brightness (inside the silhouette) and the silhouette itself, 0-1."""
    a = frame[..., 3].astype(float) / 255.0
    lum = (0.299 * frame[..., 0] + 0.587 * frame[..., 1] + 0.114 * frame[..., 2]) / 255.0
    return np.stack([lum * a, a])


def _rotate(img: np.ndarray, deg: float, c: Tuple[float, float]) -> np.ndarray:
    """Each feature plane rotated by `deg` (anticlockwise on screen) about c = (x, y), bilinear."""
    Image, _, _ = _pil()
    return np.stack([np.array(Image.fromarray(p.astype(np.float32), "F").rotate(deg, resample=Image.BILINEAR, center=c)) for p in img])


def _rot_point(x: float, y: float, deg: float, c: Tuple[float, float]) -> Tuple[float, float]:
    """Where a point goes under the same rotation as _rotate (image y points down)."""
    a = np.radians(deg); dx, dy = x - c[0], y - c[1]
    return c[0] + np.cos(a) * dx + np.sin(a) * dy, c[1] - np.sin(a) * dx + np.cos(a) * dy


def _parabola(vals: np.ndarray, i: int) -> float:
    """Sub-step position of a minimum at index i (parabola through i-1, i, i+1)."""
    if 0 < i < len(vals) - 1:
        a, b, c = vals[i - 1], vals[i], vals[i + 1]; den = a - 2 * b + c
        return i + (0.5 * (a - c) / den if den > 1e-12 else 0.0)
    return float(i)


def track_head(ref: np.ndarray, cur: np.ndarray, box: Tuple[int, int, int, int], c: Tuple[float, float],
               guess: Tuple[float, float], angles=np.arange(-24, 25, 3), reach: int = 5) -> Tuple[float, float, float]:
    """The rigid move (degrees, dx, dy) that best carries the reference frame's head (inside `box` = x0, y0, x1, y1)
    onto the current frame: rotation about c, then a shift searched within `reach` px of `guess`; sub-pixel and
    sub-step by parabola fits. Matching the whole head (ears, eyes, nose, outline) at once is far steadier than
    re-finding the ears and nose on every frame."""
    fr, fc = _features(ref), _features(cur); x0, y0, x1, y1 = box; H, W = fc.shape[1:]
    gx, gy = int(round(guess[0])), int(round(guess[1]))
    best = (np.inf, 0.0, 0.0, 0.0); by_angle = []
    for deg in angles:
        patch = _rotate(fr, float(deg), c)[:, y0:y1, x0:x1]
        grid = np.full((2 * reach + 1, 2 * reach + 1), np.inf)
        for iy, dy in enumerate(range(gy - reach, gy + reach + 1)):
            for ix, dx in enumerate(range(gx - reach, gx + reach + 1)):
                if y0 + dy < 0 or x0 + dx < 0 or y1 + dy > H or x1 + dx > W:
                    continue
                grid[iy, ix] = float(((fc[:, y0 + dy:y1 + dy, x0 + dx:x1 + dx] - patch) ** 2).sum())
        iy, ix = np.unravel_index(np.argmin(grid), grid.shape); by_angle.append(grid[iy, ix])
        if grid[iy, ix] < best[0]:
            sx = _parabola(grid[iy, :], ix) - reach + gx if np.isfinite(grid[iy, :]).all() else ix - reach + gx
            sy = _parabola(grid[:, ix], iy) - reach + gy if np.isfinite(grid[:, ix]).all() else iy - reach + gy
            best = (grid[iy, ix], float(deg), float(sx), float(sy))
    ia = int(np.argmin(by_angle)); step = float(angles[1] - angles[0])
    deg = float(angles[0] + step * _parabola(np.array(by_angle), ia))
    return deg, best[2], best[3]


def bandana(loop: np.ndarray) -> np.ndarray:
    """The dance loop (frames, h, w, 4) -> the same loop enlarged UP x, each frame wearing the bandana.

    The band is fitted ONCE, on a reference frame (the one whose ears/nose/band geometry is most typical of the loop);
    every other frame gets that same band carried by the head's rigid move (track_head), so the bandana stays put on
    the head instead of wobbling with frame-by-frame detection noise. The knot keeps its spot on the band."""
    raw = [head_geometry(f) for f in loop]
    keys = ("lx", "ly", "rx", "ry", "nose")
    med = {k: float(np.median([g[k] for g in raw])) for k in keys}
    r = int(np.argmin([sum(abs(g[k] - med[k]) for k in keys) for g in raw]))
    ref = dict(raw[r])
    # the reference band's tilt: the MEDIAN ear-line angle over the loop (one frame's ears can be tilted by its pose)
    ang = float(np.median([np.arctan2(g["ry"] - g["ly"], g["rx"] - g["lx"]) for g in raw]))
    mx, my, half = (ref["lx"] + ref["rx"]) / 2, (ref["ly"] + ref["ry"]) / 2, np.hypot(ref["rx"] - ref["lx"], ref["ry"] - ref["ly"]) / 2
    ref.update(lx=mx - half * np.cos(ang), ly=my - half * np.sin(ang), rx=mx + half * np.cos(ang), ry=my + half * np.sin(ang))
    knot_t = float(np.clip((ref["nose"] - ref["lx"]) / max(ref["rx"] - ref["lx"], 1e-6), 0.3, 0.7))
    op = loop[r][..., 3] > 128; y0 = int(np.nonzero(op.any(axis=1))[0][0])
    box = (max(0, int(ref["lx"]) - 3), max(0, y0 - 2), min(loop.shape[2], int(ref["rx"]) + 4), min(loop.shape[1], y0 + 24))
    c = ((ref["lx"] + ref["rx"]) / 2, (ref["ly"] + ref["ry"]) / 2 + 4.0)

    def head_centre(f):
        o = f[..., 3] > 128; top = int(np.nonzero(o.any(axis=1))[0][0]); ys, xs = np.nonzero(o[top:top + 22])
        return float(xs.mean()), float(ys.mean() + top)
    hc = [head_centre(f) for f in loop]
    moves = [(0.0, 0.0, 0.0) if i == r else track_head(loop[r], loop[i], box, c, (hc[i][0] - hc[r][0], hc[i][1] - hc[r][1]))
             for i in range(len(loop))]
    n = len(loop)
    # the head TURNS (3D) more than it tilts; a 2D matcher reads the turn as tilt, so only part of the fitted rotation
    # is kept (TILT_KEEP, capped at +-TILT_MAX degrees) — the band follows the head's shift fully and its tilt gently
    moves = [(float(np.clip(TILT_KEEP * m[0], -TILT_MAX, TILT_MAX)), m[1], m[2]) for m in moves]
    k = 2                                                                      # smoothing (the loop wraps)
    moves = [tuple(float(np.mean([moves[(i + d) % n][j] for d in range(-k, k + 1)])) for j in range(3)) for i in range(n)]
    geo = []
    for i, (deg, dx, dy) in enumerate(moves):
        lx, ly = _rot_point(ref["lx"], ref["ly"], deg, c); rx, ry = _rot_point(ref["rx"], ref["ry"], deg, c)
        geo.append({"lx": lx + dx, "ly": ly + dy, "rx": rx + dx, "ry": ry + dy, "thick": ref["thick"], "knot_t": knot_t,
                    "nose": ref["nose"], "cx": raw[i]["cx"], "deg": deg, "dx": dx, "dy": dy})
    for i in range(n):                                                         # sideways body speed (the knot ends' swing)
        geo[i]["vx"] = (raw[(i + 1) % n]["cx"] - raw[i - 1]["cx"]) / 2
    return np.stack([np.array(_draw(_upscale(f), g, i)) for i, (f, g) in enumerate(zip(loop, geo))])
