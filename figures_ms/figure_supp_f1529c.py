"""Supplementary Figure 2 — Vgsc F1529C in the Bouaké transcriptomes and its haplotype background.

A  RNA-seq read fractions carrying the derived allele at V402L, I1527T and F1529C in each Bouaké
   An. coluzzii sample (the V402L estimate sums its two substitutions, G>T and G>C).
B  The 402L substitution carried by phased 1527T haplotypes (Ag3 An. coluzzii, sample sets without
   a use restriction), with and without F1529C.
C  Haplotypes over +/-25 kb around codon 1529, grouped by F1529C status and 402L substitution, with
   the SNP that separates the two 402L lineages marked.

Panels B and C use only the haplotype background of F1529C; the geographic and temporal survey in
Ag3 is a separate analysis (figure_f1529c.py).
"""

from pathlib import Path

import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import pandas as pd

from f1529c_panels import RESULTS, draw_background, draw_bouake_reads, panel_haplotypes
from bouake.paths import FONT_DIR

OUT = Path(__file__).resolve().parent

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 7.5, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.edgecolor": "#4a4a4a", "svg.fonttype": "none", "figure.dpi": 150,
})


def build() -> None:
    background = pd.read_csv(RESULTS / "f1529c_haplotype_background.csv")
    fig = plt.figure(figsize=(7.2, 6.6))
    ax_reads = fig.add_axes([0.115, 0.585, 0.36, 0.345])
    ax_background = fig.add_axes([0.70, 0.725, 0.24, 0.12])
    ax_strip = fig.add_axes([0.075, 0.075, 0.05, 0.37])
    ax_heat = fig.add_axes([0.13, 0.075, 0.62, 0.37])
    draw_bouake_reads(ax=ax_reads)
    draw_background(ax=ax_background, background=background)
    panel_haplotypes(ax_heat=ax_heat, ax_strip=ax_strip, background=background)
    for letter, (x, y) in {"A": (0.02, 0.945), "B": (0.56, 0.945), "C": (0.02, 0.50)}.items():
        fig.text(x, y, letter, fontsize=13, fontweight="bold", va="bottom", ha="left")
    fig.text(0.055, 0.502, "Haplotypes around codon 1529 (±25 kb)", fontsize=8.5, fontweight="bold", va="bottom")
    fig.text(0.60, 0.947, "402L substitution on 1527T haplotypes", fontsize=8.5, fontweight="bold", va="bottom")
    fig.text(0.075, 0.947, "Bouaké RNA-seq reads", fontsize=8.5, fontweight="bold", va="bottom")
    fig.savefig(OUT / "figure_supp_f1529c.png", dpi=300)
    fig.savefig(OUT / "figure_supp_f1529c.svg")


if __name__ == "__main__":
    build()
