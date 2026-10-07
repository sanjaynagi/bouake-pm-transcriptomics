"""GSEA figure — gene-set enrichment among genes over-expressed in pirimiphos-methyl survivors.

Back-to-back bars: An. gambiae grows left, An. coluzzii right, one row per gene-set cluster, so
every cluster is tested and shown in both species. Terms that share genes (for example ABC_tran,
ABC_membrane, KEGG ABC transporters and the GO ATPase-coupled transporter term) are merged into
one cluster, whose value in each species is the best adjusted P among its member terms; this
stops one gene set being counted several times. A hollow bar is a cluster that is not
significant in that species (padj >= 0.05), drawn at its actual value.

Source: figures_ms/gsea-bouake.xlsx (hypergeometric test, over-expressed genes only).
Writes figure_gsea.{png,svg} and results/gsea/gsea_significant_terms.csv (every significant
term with its cluster, for the supplement).
"""


import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from bouake.paths import REPO, FONT_DIR

OUT = REPO / "figures_ms"
GAMBIAE, COLUZZII = "#5D69B1", "#E58606"
INK, INK_SOFT, GRID = "#1a1a1a", "#6a6a6a", "#e6e6e6"
ALPHA = 0.05
X_MAX = 21.0

for ttf in FONT_DIR.glob("FiraSans-*.ttf"):
    fm.fontManager.addfont(str(ttf))
mpl.rcParams.update({
    "font.family": "Fira Sans", "font.size": 8, "axes.labelsize": 8.5, "xtick.labelsize": 7.5,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "axes.edgecolor": "#4a4a4a",
    "svg.fonttype": "none", "figure.dpi": 150,
})

# Cluster -> (group, headline label, member annotations). Members are matched by annotation id
# across the GO, PFAM and KEGG sheets; the first member names the cluster's source term.
CLUSTERS = [
    ("Cuticle", "Chitin-binding domain", ["CBM_14"]),
    ("Transport", "ABC transporters",
     ["ABC_tran", "ABC_membrane", "ABC2_membrane_3", "aga02010", "GO:0042626"]),
    ("Transport", "Transmembrane transport",
     ["GO:0055085", "GO:0022857", "GO:0016021", "GO:0016020", "GO:1902600", "MFS_1"]),
    ("Protein homeostasis", "Ubiquitin family", ["ubiquitin"]),
    ("Protein homeostasis", "Chaperones and protein folding",
     ["HSP90", "HSP70", "GO:0051082", "GO:0034605", "HATPase_c", "aga04141"]),
    ("Protein homeostasis", "Aminopeptidases",
     ["ERAP1_C", "Peptidase_M1", "GO:0070006", "GO:0004177", "GO:0043171", "GO:0042277"]),
    ("Lipid metabolism", "Lipid metabolism", ["GO:0006629", "GO:0004806"]),
    ("Olfaction", "Odorant-binding proteins",
     ["PBP_GOBP", "GO:0005549", "GO:0007608", "GO:0042048", "GO:0035274", "GO:0035275"]),
    ("Olfaction", "DUF753 (unknown function)", ["DUF753"]),
]
SHEETS = {"GO": "go", "PFAM": "pfam", "KEGG": "kegg"}
SPECIES = {"gambiae": "gambiaeCont_gambiaePM", "coluzzii": "coluzziiCont_coluzziiPM"}
PFAM_LABEL = {
    "CBM_14": "Chitin-binding domain", "ubiquitin": "Ubiquitin family", "ABC_tran": "ABC transporter",
    "ABC_membrane": "ABC transmembrane region", "ERAP1_C": "ERAP1-like C-terminal domain",
    "Peptidase_M1": "Peptidase M1 N-terminal", "HSP90": "Hsp90 protein", "ABC2_membrane_3": "ABC-2 transporter",
    "MFS_1": "Major facilitator superfamily", "HSP70": "Hsp70 protein", "HATPase_c": "Histidine kinase-like ATPase",
    "PBP_GOBP": "Pheromone/general odorant-binding protein", "DUF753": "Domain of unknown function",
}
GO_LABEL_FIX = {"GO:0016021": "integral component of membrane", "GO:0016020": "membrane",
                "GO:0004177": "aminopeptidase activity"}


def load_all_terms() -> pd.DataFrame:
    """Every tested term in both species: annotation, database, label, padj per species."""
    wide = None
    for species, prefix in SPECIES.items():
        frames = []
        for database, suffix in SHEETS.items():
            table = pd.read_excel(OUT / "gsea-bouake.xlsx", sheet_name=f"{prefix}_{suffix}")
            table = table.drop_duplicates(subset="annotation")
            label = table["descriptions"] if "descriptions" in table else (
                table["description"] if "description" in table else table["annotation"])
            frames.append(pd.DataFrame({"annotation": table["annotation"], "database": database,
                                        "label": label, f"padj_{species}": table["padj"]}))
        part = pd.concat(frames, ignore_index=True)
        wide = part if wide is None else wide.merge(
            part[["annotation", f"padj_{species}"]], on="annotation", how="outer")
    wide["label"] = wide["annotation"].map(PFAM_LABEL).fillna(wide["label"])
    wide["label"] = wide["annotation"].map(GO_LABEL_FIX).fillna(wide["label"])
    return wide


def cluster_table(*, terms: pd.DataFrame) -> pd.DataFrame:
    """One row per cluster: best padj in each species and how many member terms are significant."""
    rows, members_long = [], []
    indexed = terms.set_index("annotation")
    for group, label, members in CLUSTERS:
        member_rows = indexed.loc[members]
        row = dict(group=group, label=label, n_members=len(members))
        for species in SPECIES:
            padj = member_rows[f"padj_{species}"]
            row[f"padj_{species}"] = padj.min()
            row[f"n_sig_{species}"] = int((padj < ALPHA).sum())
            row[f"best_{species}"] = padj.idxmin()
        rows.append(row)
        for annotation, member in member_rows.iterrows():
            members_long.append(dict(cluster=label, annotation=annotation, database=member["database"],
                                     label=member["label"], padj_gambiae=member["padj_gambiae"],
                                     padj_coluzzii=member["padj_coluzzii"]))
    return pd.DataFrame(rows), pd.DataFrame(members_long)


def draw_side(*, ax, rows: pd.DataFrame, y: np.ndarray, species: str, colour: str, side: int) -> None:
    """Bars for one species; side = -1 grows left (axis inverted by the caller), +1 right."""
    for position, (_, row) in zip(y, rows.iterrows()):
        value = -np.log10(row[f"padj_{species}"])
        significant = row[f"padj_{species}"] < ALPHA
        ax.barh(position, value, height=0.66, color=colour if significant else "none",
                edgecolor=colour, linewidth=0.8 if not significant else 0, alpha=0.95 if significant else 0.55,
                hatch=None, zorder=3)
        text = f"{row[f'padj_{species}']:.0e}".replace("e-0", "e-").replace("e+0", "e+") \
            if row[f"padj_{species}"] < 0.1 else f"n.s. ({row[f'padj_{species}']:.2f})"
        ax.text(value + 0.35, position, text, va="center", ha="left" if side > 0 else "right", fontsize=6.6,
                color=INK if significant else INK_SOFT, zorder=4)
    ax.axvline(-np.log10(ALPHA), color=INK_SOFT, lw=0.6, ls=(0, (3, 2)), zorder=2)
    ax.set_xlim(0, X_MAX)
    ax.set_xticks([0, 5, 10, 15, 20])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_yticks([])
    ax.grid(axis="x", color=GRID, lw=0.7, zorder=0)


def build() -> None:
    terms = load_all_terms()
    clusters, members = cluster_table(terms=terms)

    # Vertical layout: a header row per group, then one row per cluster.
    y, headers, cursor, previous = [], [], 0.0, None
    for group in clusters["group"]:
        if group != previous:
            cursor += 0.35 if previous is not None else 0.0
            headers.append((cursor, group))
            cursor += 0.95
            previous = group
        y.append(cursor)
        cursor += 1.0
    y = np.array(y)

    fig = plt.figure(figsize=(7.2, 0.33 * cursor + 1.25))
    grid = GridSpec(1, 3, figure=fig, width_ratios=[1.0, 0.9, 1.0], wspace=0.04,
                    left=0.03, right=0.975, top=0.915, bottom=0.17)
    ax_g, ax_mid, ax_c = (fig.add_subplot(grid[0, i]) for i in range(3))
    draw_side(ax=ax_g, rows=clusters, y=y, species="gambiae", colour=GAMBIAE, side=-1)
    draw_side(ax=ax_c, rows=clusters, y=y, species="coluzzii", colour=COLUZZII, side=1)
    ax_g.invert_xaxis()
    for ax in (ax_g, ax_c):
        ax.set_ylim(cursor - 0.3, -0.7)
    ax_mid.set_ylim(cursor - 0.3, -0.7)
    ax_mid.axis("off")

    # Centre spine: cluster label with a one-line account of what was merged.
    for position, (_, row) in zip(y, clusters.iterrows()):
        n_sig = max(row["n_sig_gambiae"], row["n_sig_coluzzii"])
        sub = f"{row['n_members']} terms merged" if row["n_members"] > 1 else "single term"
        ax_mid.text(0.5, position - 0.12, row["label"], ha="center", va="center", fontsize=8, color=INK)
        ax_mid.text(0.5, position + 0.30, sub, ha="center", va="center", fontsize=6.3, color=INK_SOFT)
    for position, group in headers:
        ax_mid.text(0.5, position, group.upper(), ha="center", va="center", fontsize=6.5,
                    color=INK_SOFT, fontweight="bold")
        ax_mid.plot([0.08, 0.92], [position + 0.28, position + 0.28], color=GRID, lw=0.8)
    ax_mid.set_xlim(0, 1)

    fig.text(0.03, 0.965, "A", fontsize=13, fontweight="bold", color=INK, va="center")
    fig.text(0.075, 0.965, "An. gambiae", fontsize=10, fontweight="bold", fontstyle="italic",
             color=GAMBIAE, va="center")
    fig.text(0.62, 0.965, "B", fontsize=13, fontweight="bold", color=INK, va="center")
    fig.text(0.665, 0.965, "An. coluzzii", fontsize=10, fontweight="bold", fontstyle="italic",
             color=COLUZZII, va="center")
    for ax in (ax_g, ax_c):
        ax.set_xlabel(r"$-\log_{10}$(adjusted $P$)")
    handles = [Line2D([], [], marker="s", ls="", color="#888", label="significant (padj < 0.05)"),
               Line2D([], [], marker="s", ls="", markerfacecolor="none", color="#888",
                      label="not significant, drawn at its value"),
               Line2D([], [], color=INK_SOFT, ls=(0, (3, 2)), lw=0.8, label="padj = 0.05")]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=7,
               bbox_to_anchor=(0.5, 0.0), handletextpad=0.4, columnspacing=1.6)

    fig.savefig(OUT / "figure_gsea.png", dpi=300)
    fig.savefig(OUT / "figure_gsea.svg")

    # Supplement: every significant term with its cluster.
    significant = terms[(terms["padj_gambiae"] < ALPHA) | (terms["padj_coluzzii"] < ALPHA)]
    significant = significant.merge(members[["annotation", "cluster"]], on="annotation", how="left")
    significant["cluster"] = significant["cluster"].fillna("unclustered")
    destination = REPO / "results" / "gsea"
    destination.mkdir(parents=True, exist_ok=True)
    significant.sort_values(["cluster", "padj_gambiae", "padj_coluzzii"]).to_csv(
        destination / "gsea_significant_terms.csv", index=False)


if __name__ == "__main__":
    build()
