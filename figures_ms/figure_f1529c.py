"""Figure 4 — Vgsc F1529C in wild An. coluzzii.

A, B  Maps of F1529C frequency in Ag3 An. coluzzii cohorts, 2012-2018 and 2019-2025.
C     F1529C frequency by year and country: absent before 2019, then rising.
D     402L substitution carried by each F1529C haplotype, by country. V402L has two independent
      origins (G>T and G>C at 2L:2,391,228); F1529C occurs on both.
E     One origin or two? Haplotype plot of 1527T haplotypes over +/-25 kb around codon 1529,
      grouped by F1529C status and 402L background.
F     The same question quantified: share of identical haplotype pairs for carrier and
      non-carrier 1527T haplotypes, holding the 402L background fixed.

Inputs: results/f1529c/ (scripts 01-04 in analysis/f1529c), Ag3 sample metadata (cached),
Natural Earth shapes in figures_ms/shp. Kept for the separate Ag3 paper on F1529C; the Bouaké
paper uses only figure_supp_f1529c.py. Cohort-level version: figure_f1529c_cohorts.py.
"""

import sys
from pathlib import Path

import malariagen_data
import matplotlib as mpl

mpl.use("Agg")  # malariagen_data registers an IPython GUI hook that pyplot would trip over
import matplotlib.font_manager as fm
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd
import shapefile
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, PowerNorm
from matplotlib.patches import Polygon

from bouake.ag3_sets import unrestricted_sample_sets
from bouake.paths import FONT_DIR, REPO

RESULTS = REPO / "results" / "f1529c"
OUT = REPO / "figures_ms"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from f1529c_panels import panel_haplotypes  # noqa: E402  (shared with the supplementary figure)

INK, INK_SOFT, GRID = "#1a1a1a", "#444444", "#e6e6e6"
LAND, LAND_EDGE = "#eceded", "#c9ccd1"
BOUAKE = (-5.032, 7.696)
PERIODS = (("2012–2018", 2012, 2018), ("2019–2025", 2019, 2025))
CMAP = LinearSegmentedColormap.from_list("f1529c", ["#fde3b8", "#E58606", "#7a3200"])
ZERO_EDGE = "#8a8a8a"
CARRIER_INK, CONTROL_GREY = "#1a1a1a", "#8a8a8a"
XLIM, YLIM = (-17.5, 6.5), (3.0, 16.5)
STUDY_COUNTRY = "Côte d'Ivoire"
COUNTRY_ORDER = ("Cote d'Ivoire", "Ghana", "Burkina Faso", "Mali", "Guinea", "Gambia, The")
COUNTRY_COLOUR = {"Burkina Faso": "#009485", "Mali": "#8f4fa8", "Ghana": "#828f1a",
                  "Cote d'Ivoire": "#c2410c", "Guinea": "#3b82a0", "Gambia, The": "#8c8c8c"}
ORIGIN_COLOUR = {"C": "#0072B2", "T": "#CC79A7", "none": "#cfcfcf"}   # Okabe-Ito blue and pink
ORIGIN_LABEL = {"C": "402L (G>C)", "T": "402L (G>T)", "none": "no 402L"}
WINDOW_BP = 25_000
# (summary key, label, group). The two KEY rows hold the 402L background constant and are the test.
CONTRASTS = (
    ("carrier[C] - carrier[T]", "F1529C carriers: 402L G>C vs G>T", "key"),
    ("control[C] - control[T]", "Non-carriers: 402L G>C vs G>T", "key"),
    ("carrier - carrier", "F1529C vs F1529C", "context"),
    ("control - control", "Non-carrier vs non-carrier", "context"),
    ("carrier - control", "F1529C vs non-carrier", "context"),
)

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 7.5, "axes.labelsize": 8, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "axes.linewidth": 0.7, "xtick.major.width": 0.7,
    "ytick.major.width": 0.7, "axes.edgecolor": "#4a4a4a", "svg.fonttype": "none", "figure.dpi": 150,
})


def tidy(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def cohort_locations() -> pd.DataFrame:
    """Mean sampling coordinates per admin1 x year, from the cached Ag3 sample metadata."""
    ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
    sample_sets = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii",))
    metadata = ag3.sample_metadata()
    metadata = metadata[(metadata["taxon"] == "coluzzii") & metadata["sample_set"].isin(sample_sets)]
    return metadata.groupby(["admin1_iso", "year"], as_index=False).agg(
        latitude=("latitude", "mean"), longitude=("longitude", "mean"))


def draw_land(*, ax) -> None:
    reader = shapefile.Reader(str(OUT / "shp" / "ne_110m_admin_0_countries.shp"))
    fields = [f[0] for f in reader.fields[1:]]
    i_name, i_sub = fields.index("NAME"), fields.index("SUBREGION")
    for record in reader.shapeRecords():
        if record.record[i_sub] not in ("Western Africa", "Middle Africa", "Northern Africa"):
            continue
        points, parts = record.shape.points, list(record.shape.parts) + [len(record.shape.points)]
        for a, b in zip(parts[:-1], parts[1:]):
            ax.add_patch(Polygon(points[a:b], closed=True, facecolor=LAND, edgecolor=LAND_EDGE,
                                 lw=0.5, zorder=1))
            if record.record[i_name] == STUDY_COUNTRY:
                ax.add_patch(Polygon(points[a:b], closed=True, facecolor="none",
                                     edgecolor=INK, lw=0.8, zorder=2))


def wilson_interval(*, count: np.ndarray, total: np.ndarray, z: float = 1.96):
    """Wilson score interval for a binomial proportion."""
    p = count / total
    centre = (p + z**2 / (2 * total)) / (1 + z**2 / total)
    half = z * np.sqrt(p * (1 - p) / total + z**2 / (4 * total**2)) / (1 + z**2 / total)
    return centre - half, centre + half


def panels_maps(*, axes, frequencies: pd.DataFrame, norm: PowerNorm) -> None:
    for ax, letter, (title, start, end) in zip(axes, "AB", PERIODS):
        draw_land(ax=ax)
        cohorts = frequencies[(frequencies["year"] >= start) & (frequencies["year"] <= end)]
        cohorts = cohorts.sort_values("frq_F1529C")    # highest drawn last
        size = np.sqrt(cohorts["n_samples"]) * 6
        absent = cohorts["count_F1529C"] == 0
        # Cohorts with no copies are open grey rings, so absence is not read as a pale low frequency.
        ax.scatter(cohorts.loc[absent, "longitude"], cohorts.loc[absent, "latitude"], s=size[absent],
                   facecolor="none", edgecolor=ZERO_EDGE, linewidth=0.6, zorder=3)
        ax.scatter(cohorts.loc[~absent, "longitude"], cohorts.loc[~absent, "latitude"], s=size[~absent],
                   c=cohorts.loc[~absent, "frq_F1529C"], cmap=CMAP, norm=norm, edgecolor="#4a4a4a",
                   linewidth=0.5, alpha=0.95, zorder=4)
        ax.plot(*BOUAKE, marker="*", ms=7, color=INK, mec="white", mew=0.5, zorder=6)
        ax.annotate("Bouaké", BOUAKE, xytext=(-14, 12), textcoords="offset points", fontsize=7, color=INK,
                    ha="right", zorder=6, path_effects=[pe.withStroke(linewidth=1.8, foreground="white")],
                    arrowprops=dict(arrowstyle="-", lw=0.6, color=INK))
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.text(0.0, 1.03, letter, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")
        ax.text(0.085, 1.045, title, transform=ax.transAxes, fontsize=8.5, fontweight="bold", va="bottom")
        ax.text(1.0, 1.045, f"{len(cohorts)} cohorts, {int(cohorts['n_samples'].sum()):,} mosquitoes",
                transform=ax.transAxes, fontsize=7, va="bottom", ha="right", color=INK_SOFT)


def map_keys(*, fig, colour_axis, size_axis, norm: PowerNorm) -> None:
    ticks = [0.01, 0.05, 0.10, 0.20, 0.35]
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=CMAP), cax=colour_axis, orientation="horizontal",
                        ticks=ticks)
    cbar.ax.set_xticklabels([f"{t:.0%}" for t in ticks])
    cbar.set_label("F1529C allele frequency (colour ramp; open grey ring = not detected)", fontsize=7, labelpad=2)
    cbar.outline.set_visible(False)
    cbar.ax.tick_params(labelsize=7, length=2)
    size_axis.axis("off")
    size_axis.text(0.0, 0.5, "Cohort size", transform=size_axis.transAxes, fontsize=7, va="center",
                   color=INK_SOFT)
    for x, size in zip((0.17, 0.42, 0.70), (50, 250, 750)):
        size_axis.scatter([x], [0.5], s=np.sqrt(size) * 6, facecolor="#f3d9b0", edgecolor="#4a4a4a",
                          linewidth=0.5, transform=size_axis.transAxes, clip_on=False)
        size_axis.text(x + 0.05, 0.5, f"{size}", transform=size_axis.transAxes, fontsize=7,
                       va="center", color=INK_SOFT)
    size_axis.text(0.86, 0.5, "mosquitoes", transform=size_axis.transAxes, fontsize=7, va="center",
                   color=INK_SOFT)


def panel_time(*, ax, frequencies: pd.DataFrame) -> None:
    """F1529C frequency by year and country (pooled cohorts), with Wilson 95% intervals."""
    cohorts = frequencies[frequencies["year"] >= 2012]
    for country in COUNTRY_ORDER:
        group = cohorts[cohorts["country"] == country]
        yearly = group.groupby("year")[["count_F1529C", "n_F1529C"]].sum()
        frequency = yearly["count_F1529C"] / yearly["n_F1529C"]
        low, high = wilson_interval(count=yearly["count_F1529C"].to_numpy(),
                                    total=yearly["n_F1529C"].to_numpy())
        colour = COUNTRY_COLOUR[country]
        for year, next_year in zip(yearly.index[:-1], yearly.index[1:]):
            if next_year - year == 1:     # join consecutive years only; gaps are not trends
                ax.plot([year, next_year], [frequency[year], frequency[next_year]], lw=1.2,
                        color=colour, alpha=0.8)
        # Intervals only where the cohort is large enough to be informative (>= 200 alleles).
        informative = (yearly["n_F1529C"] >= 200).to_numpy()
        ax.vlines(yearly.index[informative], low[informative], high[informative], color=colour,
                  lw=0.8, alpha=0.5)
        ax.scatter(yearly.index, frequency, s=np.sqrt(yearly["n_F1529C"]) * 1.1 + 4, color=colour,
                   edgecolor="white", linewidth=0.4, zorder=3, label=country.replace("Cote", "Côte").replace(", The", ""))
    ax.set_ylabel("F1529C frequency")
    ax.set_xlabel("Collection year")
    ax.set_xticks(range(2012, 2026, 2))
    ax.set_ylim(-0.012, 0.235)
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.legend(frameon=False, fontsize=6.8, loc="upper left", ncol=2, handlelength=0.8, handletextpad=0.3,
              columnspacing=0.9, borderaxespad=0.1)
    tidy(ax)


def panel_origin_by_country(*, ax, background: pd.DataFrame) -> None:
    """402L origin carried by the phased F1529C haplotypes, per country."""
    carriers = background[background["f1529c"] == 1]
    rows = [("All", carriers)] + [(c, carriers[carriers["country"] == c]) for c in COUNTRY_ORDER
                                  if (carriers["country"] == c).sum() >= 5]
    y = np.arange(len(rows))[::-1].astype(float)
    y[1:] -= 0.35                                  # gap under the pooled row
    for position, (name, subset) in zip(y, rows):
        left = 0.0
        for origin in ("C", "T", "none"):
            share = (subset["v402l_origin"] == origin).mean()
            ax.barh(position, share, left=left, color=ORIGIN_COLOUR[origin], height=0.66,
                    label=ORIGIN_LABEL[origin] if position == y[0] else None)
            if share >= 0.12:
                ax.text(left + share / 2, position, f"{share:.0%}", ha="center", va="center", fontsize=6.3,
                        color="white", fontweight="bold")
            left += share
        ax.text(1.04, position, f"n = {len(subset)}", va="center", fontsize=6.8, color=INK_SOFT)
    ax.set_yticks(y)
    ax.set_yticklabels([n.replace("Cote", "Côte").replace(", The", "") for n, _ in rows])
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("Share of phased F1529C haplotypes")
    ax.legend(frameon=False, fontsize=6.8, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3,
              handlelength=0.9, handletextpad=0.3, columnspacing=0.8, title=None)
    tidy(ax)


def panel_origin_test(*, ax, summary: pd.DataFrame) -> None:
    """Haplotype sharing around codon 1529: the key contrast holds the 402L background fixed."""
    window = summary[summary["window_bp"] == WINDOW_BP].set_index("comparison")
    y = np.array([4.0, 3.0, 1.7, 0.7, -0.3])
    for position, (key, label, group) in zip(y, CONTRASTS):
        value = window.loc[key, "share_identical"]
        carrier = key.startswith("carrier") and "control" not in key
        colour = CARRIER_INK if carrier else CONTROL_GREY
        low, high = window.loc[key, "ci_low"], window.loc[key, "ci_high"]
        ax.hlines(position, low, high, color=colour, lw=1.6, alpha=0.55, zorder=2)
        ax.scatter(value, position, s=46 if group == "key" else 30, color=colour if carrier else "white",
                   edgecolor=colour, linewidth=1.3, zorder=3)
        ax.text(high + 0.015, position, f"{value:.0%}", va="center", fontsize=7.5 if group == "key" else 7,
                color=INK, fontweight="bold" if group == "key" else "normal")
    ax.axhline(2.35, color=GRID, lw=0.8)
    ax.text(0.006, 4.72, "Same two 402L backgrounds compared", fontsize=6.8, color=INK_SOFT, va="center")
    ax.text(0.006, 2.05, "Context", fontsize=6.8, color=INK_SOFT, va="center")
    carriers_value = window.loc["carrier[C] - carrier[T]", "share_identical"]
    controls_value = window.loc["control[C] - control[T]", "share_identical"]
    ax.annotate("", xy=(carriers_value, 3.5), xytext=(controls_value, 3.5),
                arrowprops=dict(arrowstyle="<->", lw=0.8, color=INK_SOFT, shrinkA=0, shrinkB=0))
    ax.text((carriers_value + controls_value) / 2, 3.64, f"{carriers_value / controls_value:.0f}-fold",
            ha="center", va="bottom", fontsize=6.8, color=INK_SOFT)
    ax.set_yticks(y)
    ax.set_yticklabels([label for _, label, _ in CONTRASTS])
    ax.set_ylim(-0.8, 5.1)
    ax.set_xlim(0, 0.62)
    ax.set_xlabel(f"Identical 1527T haplotype pairs across ±{WINDOW_BP // 1000} kb of codon 1529")
    ax.xaxis.set_major_formatter(mticker.PercentFormatter(1.0, decimals=0))
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    tidy(ax)


def build() -> None:
    frequencies = pd.read_csv(RESULTS / "vgsc_aa_frequencies.csv")
    frequencies = frequencies[frequencies["taxon"] == "coluzzii"]
    located = frequencies.merge(cohort_locations(), on=["admin1_iso", "year"], how="left")
    located = located.dropna(subset=["latitude"])
    located = located[located["longitude"].between(*XLIM) & located["latitude"].between(*YLIM)]
    background = pd.read_csv(RESULTS / "f1529c_haplotype_background.csv")
    summary = pd.read_csv(RESULTS / "f1529c_origins_summary.csv")
    norm = PowerNorm(gamma=0.5, vmin=0, vmax=0.35)

    fig = plt.figure(figsize=(7.2, 9.4))
    # Explicit positions (figure fractions) so the equal-aspect maps leave no dead space.
    axes_maps = [fig.add_axes([0.03, 0.765, 0.46, 0.185]), fig.add_axes([0.52, 0.765, 0.46, 0.185])]
    panels_maps(axes=axes_maps, frequencies=located, norm=norm)
    colour_axis = fig.add_axes([0.12, 0.737, 0.34, 0.009])
    size_axis = fig.add_axes([0.52, 0.718, 0.46, 0.04])
    map_keys(fig=fig, colour_axis=colour_axis, size_axis=size_axis, norm=norm)

    ax_time = fig.add_axes([0.105, 0.552, 0.375, 0.135])
    ax_origin = fig.add_axes([0.665, 0.552, 0.25, 0.135])
    ax_strip = fig.add_axes([0.075, 0.255, 0.05, 0.185])
    ax_heat = fig.add_axes([0.13, 0.255, 0.62, 0.185])
    ax_test = fig.add_axes([0.31, 0.045, 0.63, 0.15])
    panel_time(ax=ax_time, frequencies=frequencies)
    panel_origin_by_country(ax=ax_origin, background=background)
    panel_haplotypes(ax_heat=ax_heat, ax_strip=ax_strip, background=background)
    panel_origin_test(ax=ax_test, summary=summary)
    for letter, (x, y) in {"C": (0.02, 0.695), "D": (0.52, 0.695), "E": (0.02, 0.475), "F": (0.02, 0.212)}.items():
        fig.text(x, y, letter, fontsize=13, fontweight="bold", va="bottom", ha="left")
    fig.text(0.055, 0.477, "Haplotypes around codon 1529 (±25 kb)", fontsize=8.5, fontweight="bold", va="bottom")
    fig.text(0.055, 0.214, "Haplotype sharing around codon 1529", fontsize=8.5, fontweight="bold", va="bottom")
    fig.savefig(OUT / "figure_f1529c.png", dpi=300)
    fig.savefig(OUT / "figure_f1529c.svg")


if __name__ == "__main__":
    build()
