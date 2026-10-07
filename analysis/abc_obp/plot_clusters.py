"""Exploratory view of the ABCH1 haplotype clustering: dendrogram with the swept cluster
marked, and its composition by country and year."""
import numpy as np, pandas as pd
import matplotlib as mpl, matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.cluster.hierarchy import dendrogram
from bouake.paths import REPO as PROJECT_ROOT, FONT_DIR

REPO = str(PROJECT_ROOT)
for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({"font.family": "Fira Sans", "font.size": 8, "svg.fonttype": "none",
                     "axes.edgecolor": "#4a4a4a", "axes.linewidth": 0.7})

Z = np.load(f"{REPO}/results/abc_obp/hapclust/linkage.npy")
h = pd.read_csv(f"{REPO}/results/abc_obp/hapclust/hap_clusters.csv")
MAIN, THRESH = 168, 0.005
SWEPT, OTHER = "#c1121f", "#b8bec9"

fig = plt.figure(figsize=(11, 7.4))
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1.0], width_ratios=[1.6, 1.0],
                      hspace=0.42, wspace=0.24)

# Dendrogram: colour only the swept clade, everything else recedes to grey.
# A link belongs to the clade when both of its children do, so propagate up from the leaves.
ax = fig.add_subplot(gs[0, :])
n = len(h)
leaf_cluster = h["cluster"].values
node_main = np.zeros(2 * n - 1, dtype=bool)
node_main[:n] = leaf_cluster == MAIN
for k in range(n - 1):
    a, b = int(Z[k, 0]), int(Z[k, 1])
    node_main[n + k] = node_main[a] and node_main[b]

dn = dendrogram(Z, no_labels=True, ax=ax,
                link_color_func=lambda k: SWEPT if node_main[k] else "#ccced2")
order = np.array(dn["leaves"])
ax.axhline(THRESH, color="#5a5a5a", lw=0.8, ls=(0, (3, 3)), zorder=5)
ax.annotate(f"cut at {THRESH} Hamming", (0.995, THRESH), xycoords=("axes fraction", "data"),
            ha="right", va="bottom", fontsize=7, color="#5a5a5a")

# Leaves sit at x = 10*i + 5, so the strip must use those coordinates
in_main = leaf_cluster[order] == MAIN
top = float(np.percentile(Z[:, 2], 99.5))
strip_y, strip_h = -0.055 * top, 0.035 * top
ax.broken_barh([(10 * i, 10) for i in np.where(in_main)[0]], (strip_y, strip_h),
               facecolors=SWEPT, edgecolors="none")
ax.set_ylim(strip_y - 0.01 * top, top)
ax.set_xlim(-20, 10 * n + 20)
ax.set_ylabel("Hamming distance")
ax.set_xticks([])
ax.spines[["top", "right", "bottom"]].set_visible(False)
ax.set_title(f"A   ABCH1 haplotype clustering — {n:,} haplotypes, West African An. coluzzii\n"
             f"     one dominant cluster (red, n={int((h['cluster']==MAIN).sum())}); "
             f"next largest = {h['cluster'].value_counts().iloc[1]}",
             loc="left", fontweight="bold", fontsize=9.5)

# Cluster size distribution
ax2 = fig.add_subplot(gs[1, 0])
sizes = h["cluster"].value_counts().sort_values(ascending=False).head(25)
cols = [SWEPT if c == MAIN else OTHER for c in sizes.index]
ax2.bar(range(len(sizes)), sizes.values, color=cols, width=0.78)
ax2.annotate(f"cluster {MAIN}\nn = {sizes.iloc[0]:,}", (0, sizes.iloc[0]), xytext=(6, -2),
             textcoords="offset points", fontsize=7.5, color=SWEPT, fontweight="bold", va="top")
ax2.set_xlabel("cluster, ranked by size (top 25)")
ax2.set_ylabel("haplotypes")
ax2.spines[["top", "right"]].set_visible(False)
ax2.set_title("B   Size distribution", loc="left", fontweight="bold", fontsize=9.5)

# Frequency of the swept cluster by country and year
ax3 = fig.add_subplot(gs[1, 1])
for country, col in [("Burkina Faso", "#5D69B1"), ("Ghana", "#E58606"), ("Mali", "#52BCA3")]:
    d = h[h["country"] == country]
    tot = d.groupby("year").size()
    n = d[d["cluster"] == MAIN].groupby("year").size().reindex(tot.index, fill_value=0)
    keep = tot >= 30
    ax3.plot(tot.index[keep], (n / tot)[keep], "o-", color=col, ms=4.5, lw=1.4, label=country)
ax3.set_ylim(0, 0.65)
ax3.set_xlabel("year")
ax3.set_ylabel(f"frequency of cluster {MAIN}")
ax3.spines[["top", "right"]].set_visible(False)
ax3.legend(frameon=False, fontsize=7.5)
ax3.set_title("C   Frequency through time", loc="left", fontweight="bold", fontsize=9.5)

fig.savefig(f"{REPO}/figures_ms/abch1_clustering_explore.png", dpi=200,
            bbox_inches="tight", facecolor="white")
print("wrote figures_ms/abch1_clustering_explore.png")
print(f"cluster {MAIN}: {(h['cluster']==MAIN).sum()} / {len(h)} haplotypes")
