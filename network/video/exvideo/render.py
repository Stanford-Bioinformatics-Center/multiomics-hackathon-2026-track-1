"""The video: fly through the figure 17 page along the walk (Node + Puppeteer), then add the song (ffmpeg)."""
from __future__ import annotations

import json
import shutil
import subprocess
from pathlib import Path
from typing import Any, Dict

from .errors import InputError, RenderError

RENDER_DIR = Path(__file__).resolve().parent.parent / "render"


def ensure_renderer() -> None:
    """Node, npm and ffmpeg must exist; install the renderer's Node packages on first use."""
    for tool in ("node", "npm", "ffmpeg"):
        if shutil.which(tool) is None:
            raise InputError(f"{tool} not found (Node.js: https://nodejs.org; ffmpeg: brew install ffmpeg)")
    if not (RENDER_DIR / "node_modules" / "puppeteer").exists():
        print("installing the renderer (puppeteer) once ...")
        res = subprocess.run(["npm", "install", "--no-audit", "--no-fund"], cwd=RENDER_DIR, capture_output=True, text=True)
        if res.returncode != 0:
            raise RenderError(f"npm install failed: {res.stderr.strip()[:500]}")


def render(spec: Dict[str, Any], audio: Path, dest: Path) -> Path:
    """Render the frames from `spec` (the timeline, overlay text and page), then mux them with the audio."""
    ensure_renderer()
    frames = dest / "frames"
    if frames.exists():
        shutil.rmtree(frames)
    frames.mkdir(parents=True)
    spec_file = dest / "render_spec.json"
    spec_file.write_text(json.dumps({**spec, "frames_dir": str(frames)}, indent=2, ensure_ascii=False), encoding="utf-8")
    res = subprocess.run(["node", str(RENDER_DIR / "render_walk.js"), str(spec_file)], cwd=RENDER_DIR)
    if res.returncode != 0:
        raise RenderError("the fly-through renderer failed (see its messages above)")
    video = dest / "music_video.mp4"
    cmd = ["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(spec["fps"]), "-i", str(frames / "%05d.jpg"), "-i", str(audio),
           "-c:v", "libx264", "-preset", "veryfast", "-crf", "20", "-pix_fmt", "yuv420p", "-c:a", "aac", "-b:a", "192k", "-shortest", str(video)]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RenderError(f"ffmpeg failed: {res.stderr.strip()[:500]}")
    shutil.rmtree(frames)
    return video
