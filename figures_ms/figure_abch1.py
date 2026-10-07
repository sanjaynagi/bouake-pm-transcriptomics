"""Figure 5 — the ABCH1 / Or39-Or38 locus: expression in Bouaké, selection in wild West
African An. coluzzii, and the architecture of the sweep.

A. ABCH1 expression across the four Bouaké contrasts, oriented so positive = higher in the
   field / resistant group. Species use the MalariaGEN Ag3 taxon colours.
B. Three selection statistics stacked over a shared genomic axis for the Ag3 An. coluzzii
   cohorts whose selection-atlas signal spans ABCH1, with a gene-model track enlarging the
   sweep summit and an inset map of the countries contributing signals.
C. Neighbour-joining tree of ABCH1 haplotypes: the swept clade is a single tight group shared
   across Burkina Faso, Ghana and Mali, so the locus carries one sweep rather than several.

Selection evidence is drawn in greyscale throughout: every cohort here is An. coluzzii, so a
taxon colour would imply a species contrast these panels do not make.

The replicated H12 apex sits ~8 kb downstream of ABCH1, between Or39 and Or38, so the gene
track labels all three genes and marks the apex rather than asserting a causal target.
"""

import os
import glob

import numpy as np
import pandas as pd
import anjl
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import shapefile
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon
import matplotlib.patheffects as pe
from bouake.paths import REPO as PROJECT_ROOT, FONT_DIR

REPO = str(PROJECT_ROOT)
SCRATCH = f"{REPO}/results/abc_obp/cache"        # gene table cached by setup.py

# MalariaGEN Ag3 TAXON_COLORS: gambiae = Vivid[1], coluzzii = Vivid[0]
GAMBIAE, COLUZZII = "#5D69B1", "#E58606"
STAT_GREY, GENE_FOCAL, GENE_OR, GENE_CTX = "#5a5e66", "#2b2d31", "#8d99ae", "#c9ccd1"
TREE_BG, TREE_SWEPT = "#d5d8dc", "#2b2d31"
INK, INK_SOFT, GRID = "#1a1a1a", "#5a5a5a", "#e6e6e6"

ABCH1 = dict(gene_id="AGAP002638", contig="2R", start=24_825_285, end=24_848_323)
XLIM = (24.30e6, 25.40e6)
ZOOM = (24.815e6, 24.872e6)
SWEPT_CLUSTER = 168

LOCUS_GENES = [("ABCH1", "AGAP002638", 1718.7), ("Or39", "AGAP002639", 6.7),
               ("Or38", "AGAP002640", 1.8)]
SIGNAL_COUNTRIES = {"Burkina Faso", "Mali", "Ghana"}
# One hue per country, shared by the statistic traces and the inset map so the map serves as
# the legend. Deliberately not orange or blue, which carry taxon meaning elsewhere in the figure.
COUNTRY_COLOUR = {"BF": "#009485", "ML": "#8f4fa8", "GH": "#828f1a"}
MAP_COLOUR = {"Burkina Faso": "#009485", "Mali": "#8f4fa8", "Ghana": "#828f1a"}
STUDY_COUNTRY = "Côte d'Ivoire"   # Natural Earth spelling; Ag3 metadata uses "Cote d'Ivoire"
BOUAKE = (-5.032, 7.696)

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({
    "font.family": "Fira Sans",
    "font.sans-serif": ["Fira Sans", "Helvetica Neue", "Arial", "DejaVu Sans"],
    "font.size": 8, "axes.labelsize": 8.5, "axes.titlesize": 9,
    "xtick.labelsize": 7.5, "ytick.labelsize": 7.5, "legend.fontsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "axes.edgecolor": "#4a4a4a", "svg.fonttype": "none", "figure.dpi": 150,
    "mathtext.fontset": "custom", "mathtext.rm": "Fira Sans",
    "mathtext.it": "Fira Sans:italic", "mathtext.bf": "Fira Sans:bold",
})

# ── Panel A ─────────────────────────────────────────────────────────────
# Kisumu_gambiaeCont / Ngousso_coluzziiCont are lab-vs-field (A vs B); their sign is flipped
# so every row reads "the field / resistant group relative to its reference".
XL = f"{REPO}/results/BouakePM_diffexp.xlsx"
CONTRASTS = [
    ("Ngousso_coluzziiCont",    "An. coluzzii", "vs Ngousso (lab)",     -1),
    ("coluzziiCont_coluzziiPM", "An. coluzzii", "PM survivor vs field", 1),
    ("Kisumu_gambiaeCont",      "An. gambiae",  "vs Kisumu (lab)",      -1),
    ("gambiaeCont_gambiaePM",   "An. gambiae",  "PM survivor vs field", 1),
]
pa = pd.DataFrame([
    dict(species=sp, label=lab,
         l2fc=sign * (r := pd.read_excel(XL, sheet_name=sh)
                      .set_index("GeneID").loc[ABCH1["gene_id"]])["log2FoldChange"],
         se=r["lfcSE"], padj=r["padj"])
    for sh, sp, lab, sign in CONTRASTS
])
pa["lo"], pa["hi"] = pa["l2fc"] - 1.96 * pa["se"], pa["l2fc"] + 1.96 * pa["se"]
stars = lambda p: "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < 0.05 else "n.s."

# ── Panel B: three statistics over the same cohorts ─────────────────────
def load_stat(pattern, key):
    out = []
    for f in sorted(glob.glob(pattern)):
        d = np.load(f, allow_pickle=True)
        x = d["x"]
        v = np.asarray(d[key], dtype=float)
        if v.ndim > 1:                       # iHS returns percentile columns
            v = np.nanmax(v, axis=1)
        n = min(len(x), len(v))
        out.append((os.path.basename(f).split("_")[0].split("-")[0], x[:n], v[:n]))
    return out

def load_pbs(path):
    """One (position, PBS) pair per cohort, from the wide table written by pbs.py."""
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c not in ("position", "n_segregating")]
    return [(c.split("-")[0], df["position"].to_numpy(), df[c].to_numpy()) for c in cols]


# G123 and iHS are computed and reported in the text but not plotted. G123 per-cohort baselines
# sit at different heights, which reads as banding rather than signal; PBS carries the same
# point as iHS but attributes the differentiation to the An. coluzzii branch specifically.
TRACKS = [
    ("H12", load_stat(f"{REPO}/results/abc_obp/h12/*.npz", "h12")),
    ("PBS", load_pbs(f"{REPO}/results/abc_obp/pbs/pbs_2R.csv")),
]

genes = pd.read_csv(f"{SCRATCH}/genes.csv")
genes = genes[genes["contig"] == ABCH1["contig"]]
zoom_genes = genes[(genes["end"] > ZOOM[0]) & (genes["start"] < ZOOM[1])]
feat = pd.read_csv(f"{REPO}/results/abc_obp/locus_features.csv")

# ── Panels C/D: ABCH1 haplotypes for both species ───────────────────────
HAP_DIR = f"{REPO}/results/abc_obp/hapclust_both"
hap = pd.read_csv(f"{HAP_DIR}/hap_meta.csv")
Zc = np.load(f"{HAP_DIR}/linkage.npy")          # complete-linkage dendrogram
Zt = np.load(f"{HAP_DIR}/njt_Z.npy")            # neighbour-joining tree
internal, leaf, edges = anjl.layout_equal_angle(Zt)

n_leaf = len(hap)
taxon = hap["taxon"].to_numpy()
SWEPT_CLUSTER = int(hap["cluster"].value_counts().idxmax())
is_swept = (hap["cluster"] == SWEPT_CLUSTER).to_numpy()

def propagate(Z, leaf_flag):
    """True for a node when every leaf beneath it satisfies leaf_flag."""
    node = np.zeros(2 * n_leaf - 1, dtype=bool)
    node[:n_leaf] = leaf_flag
    for k in range(len(Z)):
        a, b = int(Z[k, 0]), int(Z[k, 1])
        node[n_leaf + k] = node[a] and node[b]
    return node

node_swept_t = propagate(Zt, is_swept)
node_gamb_t = propagate(Zt, taxon == "gambiae")
node_colu_t = propagate(Zt, taxon == "coluzzii")

swept_counts = hap[is_swept]["country"].value_counts()
n_swept = int(is_swept.sum())
largest_gamb = int(hap[hap["taxon"] == "gambiae"]["cluster"].value_counts().iloc[0])
n_gamb = int((taxon == "gambiae").sum())
n_colu = int((taxon == "coluzzii").sum())

# ── Figure ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(7.2, 8.75))
gs = GridSpec(8, 2,
              height_ratios=[1.00, 0.46, 0.78, 0.78, 0.22, 0.56, 0.58, 1.45],
              width_ratios=[1.0, 1.0], hspace=0.0, wspace=0.10,
              left=0.235, right=0.975, top=0.960, bottom=0.045)

# ---- A: forest plot ----
axA = fig.add_subplot(gs[0, :])
ypos = np.arange(len(pa))[::-1]
for y, (_, r) in zip(ypos, pa.iterrows()):
    c = GAMBIAE if r["species"] == "An. gambiae" else COLUZZII
    axA.plot([r["lo"], r["hi"]], [y, y], color=c, lw=2.0, solid_capstyle="round", zorder=2)
    axA.plot(r["l2fc"], y, "o", ms=7.5, color=c, mec="white", mew=1.4, zorder=3)
    axA.annotate(f"{r['l2fc']:+.2f} {stars(r['padj'])}", (r["hi"], y),
                 xytext=(6, 0), textcoords="offset points", va="center", ha="left",
                 fontsize=7.5, color=INK_SOFT)
axA.axvline(0, color=INK_SOFT, lw=0.7, ls=(0, (3, 3)), zorder=1)
axA.set_yticks(ypos)
axA.set_yticklabels(pa["label"], fontsize=7.5, color=INK)
for sp, rows_, c in [("An. coluzzii", (3, 2), COLUZZII), ("An. gambiae", (1, 0), GAMBIAE)]:
    axA.annotate(sp, xy=(-0.235, np.mean(rows_)), xycoords=("axes fraction", "data"),
                 ha="left", va="center", fontsize=8, style="italic",
                 fontweight="bold", color=c, annotation_clip=False)
axA.set_xlim(-0.5, 2.15)
axA.set_ylim(-0.6, len(pa) - 0.4)
axA.set_xlabel("ABCH1 log$_2$ fold change (95% CI)", labelpad=2)
axA.spines[["top", "right", "left"]].set_visible(False)
axA.tick_params(axis="y", length=0)
axA.xaxis.grid(True, color=GRID, lw=0.6)
axA.set_axisbelow(True)
axA.set_title("A", loc="left", fontsize=11, fontweight="bold", color=INK, pad=6)

# ---- B: stacked selection statistics ----
stat_axes = []
for i, (name, series) in enumerate(TRACKS):
    ax = fig.add_subplot(gs[2 + i, :])
    ax.axvspan(*ZOOM, color="#ececea", lw=0, zorder=0)
    for code, x, v in series:
        m = (x >= XLIM[0]) & (x <= XLIM[1])
        ax.plot(x[m], v[m], color=COUNTRY_COLOUR[code], lw=0.9, alpha=0.75, zorder=3)
    ax.set_xlim(*XLIM)
    ax.set_ylim(bottom=0)
    if i:
        # Headroom so the top tick label cannot collide with the track above
        vmax = max(np.nanmax(v[(x >= XLIM[0]) & (x <= XLIM[1])]) for _, x, v in series)
        ax.set_ylim(0, vmax * 1.22)
        ax.yaxis.set_major_locator(mpl.ticker.MaxNLocator(nbins=4, steps=[1, 2, 5, 10]))
    ax.set_ylabel(name)
    ax.spines[["top", "right"]].set_visible(False)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    if i < len(TRACKS) - 1:
        ax.tick_params(labelbottom=False)
    else:
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v/1e6:.1f}"))
    if i == 0:
        ax.set_title("B", loc="left", fontsize=11, fontweight="bold", color=INK, pad=6)
    stat_axes.append(ax)
axB_top, axB_bot = stat_axes[0], stat_axes[-1]

# Inset map on the top statistic track
axM = axB_top.inset_axes([0.640, 0.0, 0.360, 1.10])
sf = shapefile.Reader(f"{REPO}/figures_ms/shp/ne_110m_admin_0_countries.shp")
fields = [f[0] for f in sf.fields[1:]]
i_name, i_sub = fields.index("NAME"), fields.index("SUBREGION")
for rec in sf.shapeRecords():
    if rec.record[i_sub] not in ("Western Africa", "Middle Africa", "Northern Africa"):
        continue
    name = rec.record[i_name]
    signal = name in SIGNAL_COUNTRIES
    pts, parts = rec.shape.points, list(rec.shape.parts) + [len(rec.shape.points)]
    for a, b in zip(parts[:-1], parts[1:]):
        axM.add_patch(Polygon(pts[a:b], closed=True,
                              facecolor=MAP_COLOUR.get(name, "#eceded"),
                              alpha=0.72 if signal else 1.0,
                              edgecolor="white" if signal else "#c9ccd1",
                              lw=0.4, zorder=3 if signal else 2))
    if name == STUDY_COUNTRY:
        for a, b in zip(parts[:-1], parts[1:]):
            axM.add_patch(Polygon(pts[a:b], closed=True, facecolor="none",
                                  edgecolor=INK, lw=0.9, zorder=4))
axM.plot(*BOUAKE, "o", ms=4.2, color=COLUZZII, mec="white", mew=0.8, zorder=6)
axM.annotate("Bouaké", BOUAKE, xytext=(-6, 0), textcoords="offset points",
             ha="right", va="center", fontsize=6.5, color=INK, zorder=6,
             path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
axM.set_xlim(-18.5, 9.5)
axM.set_ylim(2.5, 17.5)
axM.set_aspect("equal", adjustable="datalim")
axM.set_xticks([])
axM.set_yticks([])
for name, s in axM.spines.items():
    # keep the bottom edge: the panel now covers the H12 axis line, so it stands in for it
    s.set_visible(name == "bottom")
axM.set_facecolor((1, 1, 1, 0.92))
axM.set_zorder(5)



# ---- Gene-model track ----
axG = fig.add_subplot(gs[5, :])
focal_ids = {gid for _, gid, _ in LOCUS_GENES}
# Genes are split onto two rows by strand, so direction is carried by position
Y_PLUS, Y_MINUS = 0.72, 0.50
LABEL_Y = 0.38          # labels sit under both rows so the strands stay adjacent
H_CDS, H_UTR = 0.10, 0.048


def draw_gene(gene_id, colour, label=None, base=None):
    """Intron line with UTR exons thin and coding exons thick, on its strand's row."""
    g = feat[(feat["ID"] == gene_id) & (feat["type"] == "gene")].iloc[0]
    y = Y_PLUS if g["strand"] == "+" else Y_MINUS
    axG.plot([g["start"], g["end"]], [y, y], color=colour, lw=0.7, zorder=3,
             solid_capstyle="butt")
    parts = feat[feat["Parent"] == f"{gene_id}-RA"]
    for _, e in parts[parts["type"] == "exon"].iterrows():
        axG.add_patch(plt.Rectangle((e["start"], y - H_UTR), max(e["end"] - e["start"], 120),
                                    2 * H_UTR, color=colour, lw=0, zorder=5))
    for _, e in parts[parts["type"] == "CDS"].iterrows():
        axG.add_patch(plt.Rectangle((e["start"], y - H_CDS), max(e["end"] - e["start"], 120),
                                    2 * H_CDS, color=colour, lw=0, zorder=6))
    if label:
        txt = f"{label}\n{base:,.0f}" if base is not None else label
        axG.annotate(txt, (0.5 * (g["start"] + g["end"]), LABEL_Y),
                     ha="center", va="top", fontsize=7, color=INK, linespacing=1.2, zorder=7)


for _, g in zoom_genes[~zoom_genes["GeneID"].isin(focal_ids)].iterrows():
    if len(feat[(feat["ID"] == g["GeneID"]) & (feat["type"] == "gene")]):
        draw_gene(g["GeneID"], GENE_CTX)
for name, gid, base in LOCUS_GENES:
    draw_gene(gid, GENE_FOCAL if name == "ABCH1" else GENE_OR, label=name, base=base)

for y, sym in [(Y_PLUS, "+"), (Y_MINUS, "\u2212")]:
    axG.annotate(sym, xy=(-0.008, y), xycoords=("axes fraction", "data"),
                 ha="right", va="center", fontsize=8.5, color=INK_SOFT,
                 annotation_clip=False)
axG.set_xlim(*ZOOM)
axG.set_ylim(0.09, 0.86)
axG.set_yticks([])
axG.spines[["top", "right", "left"]].set_visible(False)
axG.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{v/1e6:.2f}"))
axG.set_xlabel(f"Chromosome arm {ABCH1['contig']} position (Mb)", labelpad=2)

# ---- C: dendrogram, leaves coloured by taxon ----
# Linkage is complete on Hamming distance (built in hapclust_both.py). scipy rejects
# count_sort and distance_sort together, so count_sort is used.
from scipy.cluster.hierarchy import dendrogram

axC = fig.add_subplot(gs[7, 0])
node_swept_c = propagate(Zc, is_swept)
dn = dendrogram(Zc, no_labels=True, ax=axC, count_sort=True,
                link_color_func=lambda k: TREE_SWEPT if node_swept_c[k] else "#b0b5bd")
for coll in axC.collections:
    coll.set_linewidth(0.3)
for ln in axC.get_lines():
    ln.set_linewidth(0.3)
order = np.array(dn["leaves"])

# Leaves as coloured points, replacing the separate taxon strip
leaf_x = 10 * np.arange(n_leaf) + 5
top = float(Zc[:, 2].max())
for value, colour in [("coluzzii", COLUZZII), ("gambiae", GAMBIAE)]:
    m = taxon[order] == value
    axC.scatter(leaf_x[m], np.full(m.sum(), -0.024 * top), s=3.5, c=colour,
                marker="o", linewidths=0, zorder=4, clip_on=False)
axC.set_ylim(-0.05 * top, top * 1.02)
# Fill the gridspec cell rather than forcing a square, so the row height governs
cell = axC.get_position()
axC.set_box_aspect((cell.height * fig.get_figheight()) / (cell.width * fig.get_figwidth()))
axC.set_xlim(-10 * n_leaf * 0.015, 10 * n_leaf * 1.015)
axC.set_ylabel("Hamming distance", labelpad=2)
axC.set_xticks([])
axC.spines[["top", "right", "bottom"]].set_visible(False)
axC.tick_params(axis="y", labelsize=6.5)
axC.annotate("C", xy=(0.0, 1.03), xycoords="axes fraction", ha="left", va="bottom",
             fontsize=11, fontweight="bold", color=INK, annotation_clip=False)

# ---- D: neighbour-joining tree, grey branches with taxon-coloured leaves ----
axD = fig.add_subplot(gs[7, 1])
ex, ey, eid = edges["x"].to_numpy(), edges["y"].to_numpy(), edges["id"].to_numpy()
eid = np.clip(eid, 0, 2 * n_leaf - 2)


def draw(mask, colour, lw, z):
    xs, ys = ex.copy(), ey.copy()
    xs[~mask] = np.nan
    ys[~mask] = np.nan
    axD.plot(xs, ys, color=colour, lw=lw, solid_capstyle="round", zorder=z)


draw(np.ones(len(ex), dtype=bool), "#d5d8dc", 0.3, 2)
draw(node_swept_t[eid], TREE_SWEPT, 1.0, 5)      # swept clade still marked

lx, ly, lid = leaf["x"].to_numpy(), leaf["y"].to_numpy(), leaf["id"].to_numpy().astype(int)
for value, colour in [("coluzzii", COLUZZII), ("gambiae", GAMBIAE)]:
    m = taxon[lid] == value
    axD.scatter(lx[m], ly[m], s=3.0, c=colour, marker="o", linewidths=0, zorder=4)
# Every swept haplotype is An. coluzzii, so drawing them black loses no taxon information
m = is_swept[lid]
axD.scatter(lx[m], ly[m], s=3.4, c=TREE_SWEPT, marker="o", linewidths=0, zorder=6)

finite = np.isfinite(ex) & np.isfinite(ey)
span = max(np.ptp(ex[finite]), np.ptp(ey[finite]))
cx, cy = np.mean([ex[finite].min(), ex[finite].max()]), np.mean([ey[finite].min(), ey[finite].max()])
half = span * 0.55
axD.set_xlim(cx - half, cx + half)
axD.set_ylim(cy - half, cy + half)
axD.set_aspect("equal")
axD.axis("off")
axD.annotate("D", xy=(0.0, 1.03), xycoords="axes fraction", ha="left", va="bottom",
             fontsize=11, fontweight="bold", color=INK, annotation_clip=False)
# Legend sits under the tree, outside the drawing area
axD.legend(handles=[
    Line2D([], [], color=COLUZZII, marker="o", ms=4, lw=0, label=f"An. coluzzii ({n_colu:,})"),
    Line2D([], [], color=GAMBIAE, marker="o", ms=4, lw=0, label=f"An. gambiae ({n_gamb:,})"),
    Line2D([], [], color=TREE_SWEPT, marker="o", ms=4, lw=1.6, label=f"swept clade ({n_swept:,}, all An. coluzzii)")],
    loc="upper center", bbox_to_anchor=(0.5, -0.02), frameon=False,
    fontsize=6.8, handlelength=1.4, labelspacing=0.3, borderaxespad=0.0)

for ext in ("svg", "png"):
    fig.savefig(f"{REPO}/figures_ms/figure_abch1.{ext}", dpi=400, bbox_inches="tight",
                facecolor="white")
print("wrote figures_ms/figure_abch1.svg / .png")
print({name: len(s) for name, s in TRACKS})
print(f"swept clade {n_swept}/{n_leaf}: {dict(swept_counts)}; "
      f"largest gambiae cluster {largest_gamb}/{n_gamb}")
