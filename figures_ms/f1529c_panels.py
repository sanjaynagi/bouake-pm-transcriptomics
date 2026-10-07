"""Shared panels for the Vgsc F1529C figures.

Used by figure_supp_f1529c.py (the supplementary figure in the Bouaké paper) and by
figure_f1529c.py (the Ag3 survey figure, for a separate paper).

  draw_bouake_reads   RNA-seq read fraction carrying each derived allele, per Bouaké sample
  draw_background     402L substitution carried by 1527T haplotypes with and without F1529C
  panel_haplotypes    haplotype plot over +/-25 kb of codon 1529 with the 402L lineage marker
"""


import matplotlib as mpl
import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from bouake.paths import REPO

RESULTS = REPO / "results" / "f1529c"
READS = REPO / "results" / "variantsOfInterest" / "csvs"

INK, INK_SOFT, GRID = "#1a1a1a", "#444444", "#e6e6e6"
COUNTRY_COLOUR = {"Burkina Faso": "#009485", "Mali": "#8f4fa8", "Ghana": "#828f1a",
                  "Cote d'Ivoire": "#c2410c", "Guinea": "#3b82a0", "Gambia, The": "#8c8c8c"}
ORIGIN_COLOUR = {"C": "#0072B2", "T": "#CC79A7", "none": "#cfcfcf"}   # Okabe-Ito blue and pink
ORIGIN_LABEL = {"C": "402L (G>C)", "T": "402L (G>T)", "none": "no 402L"}
CHANGE_COLOUR = {"V402L": "#7a7a7a", "I1527T": "#5D69B1", "F1529C": "#E58606"}
# derived nucleotides per site (V402L is two substitutions, G>T and G>C, both encoding leucine)
DERIVED = {"V402L": ("T", "C"), "I1527T": ("C",), "F1529C": ("G",)}


def tidy(ax) -> None:
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def allele_balance(*, change: str) -> pd.DataFrame:
    """Derived-read fraction, depth and 95% Wilson interval per sample at one site."""
    table = pd.read_csv(READS / f"{change}_alleleBalance.csv")
    derived = table[list(DERIVED[change])].sum(axis=1)
    depth = table[["A", "C", "G", "T"]].sum(axis=1)
    fraction = derived / depth
    z = 1.96
    centre = (fraction + z**2 / (2 * depth)) / (1 + z**2 / depth)
    half = z * np.sqrt(fraction * (1 - fraction) / depth + z**2 / (4 * depth**2)) / (1 + z**2 / depth)
    return pd.DataFrame({"sample": table["sample"], "change": change, "fraction": fraction,
                         "depth": depth, "low": centre - half, "high": centre + half})


def draw_bouake_reads(*, ax) -> None:
    """Per-sample read fractions at V402L, I1527T and F1529C in the Bouaké An. coluzzii samples."""
    balance = pd.concat([allele_balance(change=c) for c in DERIVED], ignore_index=True)
    samples = [f"coluzziiCont{i}" for i in range(1, 5)] + [f"coluzziiPM{i}" for i in range(1, 5)]
    labels = [f"Control {i}" for i in range(1, 5)] + [f"Survivor {i}" for i in range(1, 5)]
    offsets = {"V402L": 0.24, "I1527T": 0.0, "F1529C": -0.24}
    for y, sample in enumerate(samples):
        for change, offset in offsets.items():
            row = balance[(balance["sample"] == sample) & (balance["change"] == change)].iloc[0]
            ax.plot([row["low"], row["high"]], [y + offset] * 2, color=CHANGE_COLOUR[change], lw=1.3,
                    solid_capstyle="round", alpha=0.5, zorder=2)
            ax.scatter(row["fraction"], y + offset, s=22, color=CHANGE_COLOUR[change], edgecolor="white",
                       linewidth=0.5, zorder=3, label=change if y == 0 else None)
    ax.axhline(3.5, color=GRID, lw=0.9)
    ax.set_yticks(range(len(samples)))
    ax.set_yticklabels(labels)
    ax.set_ylim(len(samples) - 0.45, -0.6)
    ax.set_xlim(-0.02, 1.0)
    ax.set_xlabel("Fraction of RNA-seq reads with the derived allele (95% CI)")
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)
    ax.legend(frameon=False, fontsize=6.8, loc="upper right", handletextpad=0.2)
    tidy(ax)


def draw_background(*, ax, background: pd.DataFrame) -> None:
    """402L substitution on 1527T haplotypes, split by whether they carry F1529C."""
    groups = (("F1529C", background[background["f1529c"] == 1]),
              ("No F1529C", background[background["f1529c"] == 0]))
    y = np.array([1.0, 0.0])
    for position, (label, subset) in zip(y, groups):
        left = 0.0
        for origin in ("C", "T", "none"):
            share = (subset["v402l_origin"] == origin).mean()
            ax.barh(position, share, left=left, color=ORIGIN_COLOUR[origin], height=0.6,
                    label=ORIGIN_LABEL[origin] if position == y[0] else None)
            if share >= 0.1:
                ax.text(left + share / 2, position, f"{share:.0%}", ha="center", va="center", fontsize=6.8,
                        color="white", fontweight="bold")
            left += share
        ax.text(1.03, position, f"n = {len(subset):,}", va="center", fontsize=6.8, color=INK_SOFT)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{label}\nhaplotypes" for label, _ in groups])
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("Share of phased 1527T haplotypes")
    ax.legend(frameon=False, fontsize=6.8, loc="lower center", bbox_to_anchor=(0.45, 1.0), ncol=3,
              handlelength=0.9, handletextpad=0.3, columnspacing=0.8)
    tidy(ax)


def panel_haplotypes(*, ax_heat, ax_strip, background: pd.DataFrame) -> None:
    """Haplotype plot of 1527T haplotypes over +/-25 kb of codon 1529.

    Rows are haplotypes, grouped as 1529C carriers on each 402L background and non-carriers on each,
    and ordered by similarity within a group; columns are phased SNPs (1527 and 1529 excluded).
    A dark cell is a difference from the commonest haplotype. If F1529C arose once, both carrier
    groups share one haplotype around codon 1529, as the two non-carrier groups, which differ
    in 402L background, do not. The five F1529C haplotypes without 402L are omitted.
    """
    saved = np.load(RESULTS / "f1529c_origins_haplotypes.npz")
    meta = background.iloc[saved["subset"]].reset_index(drop=True)
    alleles, positions = saved["alleles"], saved["positions"]
    # The 1527T lineage is a recent sweep, so most SNPs are invariant: plot the variable ones
    # (at least two copies of the minor allele among the plotted haplotypes).
    minor_count = np.minimum(alleles.sum(axis=0), (1 - alleles).sum(axis=0))
    keep = minor_count >= 2
    alleles, positions = alleles[:, keep], positions[keep]
    # Polarise to the commonest haplotype: a dark cell is a difference from it, so the shared
    # sweep haplotype is blank and each haplotype's departures from it stand out.
    alleles = (alleles != (alleles.mean(axis=0) > 0.5)).astype(np.int8)

    carrier = (meta["f1529c"] == 1).to_numpy()
    origin = meta["v402l_origin"].to_numpy()
    groups = [
        ("F1529C on 402L G>C", carrier & (origin == "C")),
        ("F1529C on 402L G>T", carrier & (origin == "T")),
        ("No F1529C, 402L G>C", ~carrier & (origin == "C")),
        ("No F1529C, 402L G>T", ~carrier & (origin == "T")),
    ]
    order, cuts = [], []
    for label, mask in groups:
        rows = np.where(mask)[0]
        if len(rows) > 2:
            rows = rows[hierarchy.leaves_list(hierarchy.linkage(alleles[rows], method="average",
                                                                 metric="hamming"))]
        cuts.append((label, len(order), len(order) + len(rows)))
        order.extend(rows)
    order = np.array(order)
    gap = 6                                      # blank rows between groups
    matrix = np.full((len(order) + gap * (len(groups) - 1), alleles.shape[1]), np.nan)
    y_mid, cursor = [], 0
    for label, start, stop in cuts:
        matrix[cursor:cursor + stop - start] = alleles[order[start:stop]]
        y_mid.append((label, cursor, cursor + stop - start))
        cursor += stop - start + gap
    cmap = mpl.colors.ListedColormap(["#f4f4f4", "#2b2d31"])
    cmap.set_bad("white")
    ax_heat.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest", vmin=0, vmax=1)
    marker = pd.read_csv(RESULTS / "f1529c_lineage_marker.csv").set_index("group")
    marker_column = int(np.searchsorted(positions, marker["position"].iloc[0]))
    share = marker["share_with_lineage_allele"]
    carriers_share = (share["F1529C carriers, 402L G>T"] * marker.loc["F1529C carriers, 402L G>T", "n_haplotypes"]
                      + share["F1529C carriers, 402L G>C"] * marker.loc["F1529C carriers, 402L G>C", "n_haplotypes"]
                      ) / (marker.loc["F1529C carriers, 402L G>T", "n_haplotypes"]
                           + marker.loc["F1529C carriers, 402L G>C", "n_haplotypes"])
    ax_heat.annotate(
        f"402L G>T lineage allele at {marker['offset_bp'].iloc[0] / 1000:+.1f} kb: carried by "
        f"{share['non-carriers, 402L G>T']:.0%} of non-carriers on G>T, {share['non-carriers, 402L G>C']:.0%} on G>C, "
        f"{carriers_share:.0%} of F1529C haplotypes",
        xy=(marker_column, 1.0), xytext=(marker_column, 1.075), xycoords=ax_heat.get_xaxis_transform(),
        textcoords=ax_heat.get_xaxis_transform(), ha="left", va="bottom", fontsize=6.6, color=INK,
        arrowprops=dict(arrowstyle="-|>", lw=0.7, color=INK, shrinkA=0, shrinkB=0, mutation_scale=7))
    focal = np.searchsorted(positions, 2429623)
    ax_heat.axvline(focal - 0.5, color="#E58606", lw=0.9)
    ticks = [np.searchsorted(positions, 2429623 + kb * 1000) for kb in (-20, -10, 0, 10, 20)]
    ax_heat.set_xticks(ticks)
    ax_heat.set_xticklabels(["−20", "−10", "codon 1529", "+10", "+20 kb"])
    ax_heat.set_yticks([])
    for spine in ("top", "right", "left"):
        ax_heat.spines[spine].set_visible(False)
    for label, start, stop in y_mid:
        ax_heat.text(1.01, (start + stop) / 2, f"{label}\n(n = {stop - start})", transform=ax_heat.get_yaxis_transform(),
                     va="center", fontsize=6.6, color=INK)

    # Strip: 402L background and country of each row, in the same row order.
    strip = np.full((matrix.shape[0], 2, 3), 1.0)
    cursor = 0
    for label, start, stop in cuts:
        rows = order[start:stop]
        for i, row in enumerate(rows):
            strip[cursor + i, 0] = mpl.colors.to_rgb(ORIGIN_COLOUR[origin[row]])
            strip[cursor + i, 1] = mpl.colors.to_rgb(COUNTRY_COLOUR.get(meta["country"].iloc[row], "#bdbdbd"))
        cursor += stop - start + gap
    ax_strip.imshow(strip, aspect="auto", interpolation="nearest")
    ax_strip.set_xticks([0, 1])
    ax_strip.set_xticklabels(["402L", "country"], rotation=90, fontsize=6.3)
    ax_strip.xaxis.tick_top()
    ax_strip.tick_params(length=0)
    ax_strip.set_yticks([])
    for spine in ax_strip.spines.values():
        spine.set_visible(False)
