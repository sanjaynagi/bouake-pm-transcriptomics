"""Supplementary figure — Vgsc allele frequencies in every An. coluzzii cohort with F1529C.

One row per admin1 x year cohort (>= 10 samples) carrying F1529C, ordered by its frequency, with
V402L (G>T plus G>C), I1527T and F1529C frequencies. F1529C is never above I1527T.

Input: results/f1529c/vgsc_aa_frequencies.csv (script 01 in analysis/f1529c).
"""


import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from bouake.paths import REPO, FONT_DIR

RESULTS = REPO / "results" / "f1529c"
OUT = REPO / "figures_ms"
INK, GRID = "#1a1a1a", "#e6e6e6"
CHANGE_COLOUR = {"V402L": "#a9a9a9", "I1527T": "#5D69B1", "F1529C": "#E58606"}

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5, "axes.linewidth": 0.7, "xtick.major.width": 0.7,
    "ytick.major.width": 0.7, "axes.edgecolor": "#4a4a4a", "svg.fonttype": "none", "figure.dpi": 150,
})


def tidy(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def panel_a(*, ax, frequencies: pd.DataFrame) -> None:
    """Cohorts with any F1529C, as paired dots per change."""
    cohorts = frequencies[(frequencies["taxon"] == "coluzzii") & (frequencies["count_F1529C"] > 0)]
    cohorts = cohorts.sort_values("frq_F1529C").reset_index(drop=True)
    y = np.arange(len(cohorts))
    for change, marker in (("V402L", "s"), ("I1527T", "o"), ("F1529C", "o")):
        ax.scatter(cohorts[f"frq_{change}"], y, s=22 if change == "F1529C" else 16,
                   color=CHANGE_COLOUR[change], marker=marker, zorder=3, label=change,
                   edgecolor="white", linewidth=0.5)
    ax.hlines(y, cohorts["frq_F1529C"], cohorts["frq_I1527T"], color=GRID, lw=1.2, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.admin1_iso}  {int(r.year)}" for r in cohorts.itertuples()], fontsize=6.3)
    ax.set_xlim(-0.02, 0.72)
    ax.set_xlabel("Allele frequency")
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.legend(frameon=False, loc="lower right", fontsize=7.5, handletextpad=0.2)
    tidy(ax)


def build() -> None:
    frequencies = pd.read_csv(RESULTS / "vgsc_aa_frequencies.csv")
    fig, ax = plt.subplots(figsize=(4.6, 6.4))
    panel_a(ax=ax, frequencies=frequencies)
    fig.tight_layout()
    fig.savefig(OUT / "figure_f1529c_cohorts.png", dpi=300)
    fig.savefig(OUT / "figure_f1529c_cohorts.svg")


if __name__ == "__main__":
    build()
