#!/usr/bin/env bash
# =====================================================================================================
# video/setup.sh — one-time set-up of the music-video tools (about 2 minutes)
# =====================================================================================================
# Creates network/video/.venv with faster-whisper (syncs the lyrics to the song; the model, ~0.5 GB, downloads on
# first use), installs the renderer's Node packages (Puppeteer + its Chrome), and checks ffmpeg and Claude Code.
# make_music_video.py switches to .venv by itself, so afterwards just run: python3 network/video/make_music_video.py
# =====================================================================================================
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
command -v ffmpeg >/dev/null || { echo "ffmpeg missing: brew install ffmpeg"; exit 1; }
command -v node >/dev/null || { echo "Node.js missing: https://nodejs.org"; exit 1; }
command -v claude >/dev/null || echo "note: Claude Code ('claude') not found; use --backend api with ANTHROPIC_API_KEY"
[ -x "$HERE/.venv/bin/python" ] || python3 -m venv "$HERE/.venv"
"$HERE/.venv/bin/pip" install -q --upgrade pip
"$HERE/.venv/bin/pip" install -q "faster-whisper==1.2.1" "pillow>=10"
(cd "$HERE/render" && npm install --no-audit --no-fund)
echo "ready: python3 $HERE/make_music_video.py"
