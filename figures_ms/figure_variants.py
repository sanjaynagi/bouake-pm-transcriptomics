"""Variant analysis figure (manuscript Figure 3): ancestry, karyotype and variants of interest.

Left, in two grey panels: ancestry (AIM proportions) and karyotype (2La outer ring, 2Rb inner ring)
for each treatment group, as donuts. Right: heatmap of variant allele frequencies at known
insecticide-resistance variants, per treatment group. Layout follows the original figure
(figure4.png, assembled by hand in Inkscape); donuts, heatmap and annotations are redrawn here as
vectors from the data.

Inputs
  results/variantsOfInterest/csvs/mean_<variant>_alleleBalance.csv   per-treatment mean read frequencies
  figures_ms/2La-freq.csv, 2Rb-freq.csv                             mean karyotype frequencies
  AIM proportions                                                   Supplementary Table (below)

A heatmap cell is the proportion of reads carrying the derived allele, averaged over the samples of a
group; the two V402L substitutions (G>T, G>C) are summed. An asterisk marks a non-zero cell with a
mean read depth below 10.
"""


import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import FancyBboxPatch, Polygon
from bouake.paths import REPO, FONT_DIR

OUT = REPO / "figures_ms"
CSVS = REPO / "results" / "variantsOfInterest" / "csvs"

WIDTH_PX, HEIGHT_PX, DPI = 3453, 1843, 300       # canvas of the original figure
PT = 72 / DPI                                     # points per canvas pixel
PANELS = {"ancestry": (40, 895), "karyotype": (925, 1830)}   # top, bottom of the grey panels (px)
PANEL_LEFT, PANEL_RIGHT = 37, 1567
GAMBIAE, COLUZZII = "#5D69B1", "#E58606"          # MalariaGEN Ag3 taxon colours
KARYO = {"2La_standard": "#52bca3", "2La_inverted": "#99c945",
         "2Rb_standard": "#cc61b0", "2Rb_inverted": "#24796c"}
PANEL, INK, SOFT = "#f4f4f4", "#222222", "#666666"
TREATMENTS = ("Kisumu", "Ngousso", "coluzziiCont", "coluzziiPM", "gambiaeCont", "gambiaePM")
# Donut grid: gambiae row on top (lab strain, control, PM), coluzzii row below.
DONUT_ORDER = (("Kisumu", "gambiaeCont", "gambiaePM"), ("Ngousso", "coluzziiCont", "coluzziiPM"))
# Supplementary Table: ancestry informative marker proportions (gambiae, coluzzii)
AIMS = {"Kisumu": (0.639, 0.361), "Ngousso": (0.062, 0.938), "coluzziiPM": (0.113, 0.887),
        "coluzziiCont": (0.124, 0.876), "gambiaeCont": (0.990, 0.010), "gambiaePM": (0.979, 0.021)}

# Heatmap rows, in the original order: (gene label, variant).
ROWS = (
    ("ACE1", "G280S"), ("VGSC", "R254K"), ("VGSC", "V402L"), ("VGSC", "D466H"), ("VGSC", "T791M"),
    ("VGSC", "L995S"), ("VGSC", "L995F"), ("VGSC", "I1527T"), ("VGSC", "F1529L"), ("VGSC", "F1529C"),
    ("VGSC", "N1570Y"), ("VGSC", "A1746S"), ("VGSC", "V1853I"), ("VGSC", "I1868T"), ("VGSC", "P1874S"),
    ("VGSC", "P1874L"), ("VGSC", "F1920S"), ("VGSC", "I1940T"), ("Rdl", "A296G"), ("Rdl", "A296S"),
    ("GSTe2", "I114T"), ("GSTe2", "L119V"), ("GSTe2", "F120L"), ("COEAE1F", "E477V"),
    ("COEAE2F", "S357N"), ("CYP4J5", "L43F"), ("CYP6P4", "I236M"), ("CYP6AA1", "D155N"),
    ("CYP6P3", "E205D"), ("CYP9K1", "N224I"), ("RDGA", "S870P"),
)
LOW_DEPTH = 10

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({"font.family": "Fira Sans", "svg.fonttype": "none", "figure.dpi": 150})


def heatmap_table() -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Mean derived-read frequency and mean depth per variant and treatment, plus row labels."""
    values, depths, labels = {}, {}, []
    for gene, variant in ROWS:
        table = pd.read_csv(CSVS / f"mean_{variant}_alleleBalance.csv").set_index("treatment")
        proportion_columns = [c for c in table.columns if c.startswith("proportion")]
        values[variant] = table.loc[list(TREATMENTS), proportion_columns].sum(axis=1)
        depths[variant] = table.loc[list(TREATMENTS), "cov_mean"]
        position = f"{table['chrom'].iloc[0]}:{int(table['pos'].iloc[0])}"
        labels.append(f"{gene} | {variant} | {position}")
    order = [variant for _, variant in ROWS]
    return pd.DataFrame(values).T.loc[order], pd.DataFrame(depths).T.loc[order], labels


def annotation(value: float, depth: float) -> str:
    """'.77' style text; blank for zero, '1' for one, an asterisk when read depth is low."""
    if value < 0.005:
        return ""
    text = "1" if value >= 0.995 else f"{value:.2f}".lstrip("0")
    return text + ("*" if depth < LOW_DEPTH else "")


def sector(*, cx, cy, r_outer, r_inner, start, stop, colour, ax) -> None:
    """Annular sector from `start` to `stop` fraction of a turn, clockwise from the top (pixel y down)."""
    if stop - start <= 1e-9:
        return
    phi = np.linspace(start, stop, max(8, int(400 * (stop - start)))) * 2 * np.pi
    outer = np.column_stack([cx + r_outer * np.sin(phi), cy - r_outer * np.cos(phi)])
    inner = np.column_stack([cx + r_inner * np.sin(phi[::-1]), cy - r_inner * np.cos(phi[::-1])])
    ax.add_patch(Polygon(np.vstack([outer, inner]), closed=True, facecolor=colour, edgecolor="none", lw=0))


def donut(*, ax, cx, cy, r_outer, r_inner, fractions, colours) -> None:
    start = 0.0
    for fraction, colour in zip(fractions, colours):
        sector(cx=cx, cy=cy, r_outer=r_outer, r_inner=r_inner, start=start, stop=start + fraction,
               colour=colour, ax=ax)
        start += fraction


def arc_label(*, ax, fig, text, cx, cy, radius, start_deg, size_pt) -> None:
    """Italic text laid along a circle, starting `start_deg` clockwise from the top."""
    renderer = fig.canvas.get_renderer()
    probe = [ax.text(0, 0, text[:i + 1], fontsize=size_pt, fontstyle="italic") for i in range(len(text))]
    widths = [p.get_window_extent(renderer).width * (WIDTH_PX / fig.bbox.width) for p in probe]
    for p in probe:
        p.remove()
    previous = 0.0
    for char, width in zip(text, widths):
        centre = (previous + width) / 2
        previous = width
        phi = start_deg + np.degrees(centre / radius)
        x, y = cx + radius * np.sin(np.radians(phi)), cy - radius * np.cos(np.radians(phi))
        ax.text(x, y, char, rotation=-phi, ha="center", va="baseline", rotation_mode="anchor",
                fontsize=size_pt, fontstyle="italic", color=SOFT)


def legend_ring(*, ax, x, y, colour, label, size_pt) -> None:
    donut(ax=ax, cx=x, cy=y, r_outer=42, r_inner=21, fractions=[1.0], colours=[colour])
    ax.text(x, y + 56, label, ha="center", va="top", fontsize=size_pt, color=SOFT)


def check_text_inside_panels(*, fig, canvas) -> None:
    """Every text artist whose anchor lies in a grey panel must sit fully inside it."""
    renderer = fig.canvas.get_renderer()
    to_canvas = canvas.transData.inverted()
    problems = []
    for text in canvas.texts:
        x0, y0, x1, y1 = to_canvas.transform(text.get_window_extent(renderer)).ravel()[[0, 1, 2, 3]]
        left, right = sorted((x0, x1))
        top, bottom = sorted((y0, y1))
        for name, (panel_top, panel_bottom) in PANELS.items():
            anchor_x, anchor_y = text.get_position()
            if PANEL_LEFT < anchor_x < PANEL_RIGHT and panel_top < anchor_y < panel_bottom:
                margin = 12
                if left < PANEL_LEFT + margin or right > PANEL_RIGHT - margin or top < panel_top + margin \
                        or bottom > panel_bottom - margin:
                    problems.append((name, text.get_text(), round(left), round(top), round(right), round(bottom)))
    if problems:
        raise ValueError(f"text outside its grey panel: {problems}")


def build() -> None:
    fig = plt.figure(figsize=(WIDTH_PX / DPI, HEIGHT_PX / DPI), dpi=DPI)
    canvas = fig.add_axes([0, 0, 1, 1])
    canvas.set_xlim(0, WIDTH_PX)
    canvas.set_ylim(HEIGHT_PX, 0)
    canvas.axis("off")

    # Grey panels
    for top, bottom in PANELS.values():
        canvas.add_patch(FancyBboxPatch((PANEL_LEFT, top), PANEL_RIGHT - PANEL_LEFT, bottom - top,
                                        boxstyle="round,pad=0,rounding_size=75", facecolor=PANEL, edgecolor="none"))

    columns = (433, 803, 1175)
    # ---- Ancestry
    rows_y = (275, 597)
    for row, (y, names) in enumerate(zip(rows_y, DONUT_ORDER)):
        for x, name in zip(columns, names):
            gambiae, coluzzii = AIMS[name]
            donut(ax=canvas, cx=x, cy=y, r_outer=137, r_inner=83, fractions=[gambiae, coluzzii],
                  colours=[GAMBIAE, COLUZZII])
    canvas.text(112, sum(PANELS["ancestry"]) / 2, "Ancestry", rotation=90, ha="center", va="center", fontsize=60 * PT, color=INK)
    legend_ring(ax=canvas, x=615, y=780, colour=GAMBIAE, label="gambiae alleles", size_pt=44 * PT)
    legend_ring(ax=canvas, x=990, y=780, colour=COLUZZII, label="coluzzii alleles", size_pt=44 * PT)

    # ---- Karyotype
    frequencies = {inv: pd.read_csv(OUT / f"{inv}-freq.csv", sep="\t", index_col=0).set_index("treatment")
                   for inv in ("2La", "2Rb")}
    rows_y_karyo = (1150, 1500)
    for y, names in zip(rows_y_karyo, DONUT_ORDER):
        for x, name in zip(columns, names):
            donut(ax=canvas, cx=x, cy=y, r_outer=136, r_inner=84,
                  fractions=[frequencies["2La"].loc[name, "standard"], frequencies["2La"].loc[name, "inverted"]],
                  colours=[KARYO["2La_standard"], KARYO["2La_inverted"]])
            donut(ax=canvas, cx=x, cy=y, r_outer=84, r_inner=50,
                  fractions=[frequencies["2Rb"].loc[name, "standard"], frequencies["2Rb"].loc[name, "inverted"]],
                  colours=[KARYO["2Rb_standard"], KARYO["2Rb_inverted"]])
    canvas.text(112, sum(PANELS["karyotype"]) / 2, "Karyotype", rotation=90, ha="center", va="center", fontsize=60 * PT, color=INK)
    for x, key, label in ((268, "2La_standard", "2La standard"), (592, "2La_inverted", "2La inverted"),
                          (953, "2Rb_standard", "2Rb standard"), (1296, "2Rb_inverted", "2Rb inverted")):
        legend_ring(ax=canvas, x=x, y=1702, colour=KARYO[key], label=label, size_pt=44 * PT)

    # Panel letters: inside the top-left corner of each grey panel, and above the heatmap labels
    for letter, (x, y) in {"A": (95, PANELS["ancestry"][0] + 62), "B": (95, PANELS["karyotype"][0] + 62),
                           "C": (1612, 46)}.items():
        canvas.text(x, y, letter, ha="center", va="center", fontsize=74 * PT, fontweight="bold", color="black")

    # Strain names on an arc, and species labels on the right of each panel
    for ys in (rows_y, rows_y_karyo):
        for y, strain in zip(ys, ("Kisumu", "Ngousso")):
            arc_label(ax=canvas, fig=fig, text=strain, cx=columns[0], cy=y, radius=178, start_deg=-84,
                      size_pt=48 * PT)
        for y, species in zip(ys, ("An. gambiae", "An. coluzzii")):
            canvas.text(1448, y, species, rotation=-90, ha="center", va="center", fontsize=48 * PT,
                        fontstyle="italic", color=SOFT)

    # ---- Heatmap
    values, depths, labels = heatmap_table()
    notes = pd.DataFrame([[annotation(values.iloc[i, j], depths.iloc[i, j]) for j in range(values.shape[1])]
                          for i in range(values.shape[0])])
    left, right, top, bottom = 2162, 3418, 83, 1640
    heat = fig.add_axes([left / WIDTH_PX, 1 - bottom / HEIGHT_PX, (right - left) / WIDTH_PX,
                         (bottom - top) / HEIGHT_PX])
    sns.heatmap(values.to_numpy(), ax=heat, cmap="rocket_r", vmin=0, vmax=1, cbar=False,
                annot=notes.to_numpy(), fmt="", annot_kws={"fontsize": 43 * PT}, linewidths=3 * PT,
                linecolor="white", xticklabels=False, yticklabels=False)
    heat.set_yticks(np.arange(len(labels)) + 0.5)
    heat.set_yticklabels(labels, fontsize=40 * PT, color="black")
    heat.set_xticks(np.arange(6) + 0.5)
    heat.set_xticklabels(["Kisumu", "Ngousso", "Control", "PM", "Control", "PM"], fontsize=44 * PT, color="black")
    heat.tick_params(axis="both", length=6, width=1.0, pad=8)
    heat.tick_params(axis="x", pad=4)
    for spine in heat.spines.values():
        spine.set_visible(False)

    # Species labels under the columns
    cell = (right - left) / 6
    for first, colour, name in ((2, COLUZZII, "An. coluzzii"), (4, GAMBIAE, "An. gambiae")):
        x0 = left + first * cell + 12
        canvas.add_patch(FancyBboxPatch((x0, 1758), 2 * cell - 24, 70, boxstyle="round,pad=0,rounding_size=22",
                                        facecolor=colour, edgecolor="none"))
        canvas.text(x0 + cell - 12, 1793, name, ha="center", va="center", fontsize=54 * PT, fontstyle="italic",
                    color="white")

    check_text_inside_panels(fig=fig, canvas=canvas)
    fig.savefig(OUT / "figure_variants.png", dpi=DPI)
    fig.savefig(OUT / "figure_variants.svg")


if __name__ == "__main__":
    build()
