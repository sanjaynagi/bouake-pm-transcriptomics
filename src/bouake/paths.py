"""Project locations, resolved from this file so that no script hardcodes a user's home directory."""

from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
RESULTS = REPO / "results"
FIGURES = REPO / "figures_ms"
FONT_DIR = Path.home() / "Library" / "Fonts"
