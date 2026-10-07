"""Shared setup for Phase B follow-up analyses: DE tables, family definitions, gene coordinates, sweep overlap."""

import os
import numpy as np
import pandas as pd

from bouake.paths import REPO

os.chdir(REPO)
# Gene annotations, gene coordinates and per-gene sweep overlap are cached here (regenerated if absent).
SCRATCH = "results/abc_obp/cache"
os.makedirs(SCRATCH, exist_ok=True)

RNG_SEED = 42
N_PERM = 10_000

# ── Bouaké DE ───────────────────────────────────────────────────────────
DE_PATH = "results/BouakePM_diffexp.xlsx"
gam = pd.read_excel(DE_PATH, sheet_name="gambiaeCont_gambiaePM")
col = pd.read_excel(DE_PATH, sheet_name="coluzziiCont_coluzziiPM")

# ── Family definitions (PFAM ∪ GO ∪ name regex), matching notebook 06 ────
import anoexpress as ae

_ANN_CACHE = f"{SCRATCH}/ann_by_gene.csv"
if os.path.exists(_ANN_CACHE):
    ann_by_gene = pd.read_csv(_ANN_CACHE)
else:
    ann = ae.load_annotations()
    ann_by_gene = (
        ann.groupby("gene_id")
        .agg(
            pfam=("pfamid", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v))))),
            go=("GO_terms", lambda x: ";".join(sorted(set(
                t.strip() for v in x if pd.notna(v) for t in str(v).split(";") if t.strip())))),
        )
        .reset_index()
    )
    ann_by_gene.to_csv(_ANN_CACHE, index=False)

FAMILIES = {
    # ABC family = the ABC nucleotide-binding domain (PF00005) or an Abc* gene name. An earlier
    # definition also took GO:0042626 and four transmembrane-domain PFAMs, which pulled in P-type and
    # V-type ATPases and other non-ABC ATP-binding proteins (about half of the set).
    "ABC":  dict(pfam=["PF00005"], go=[], name_regex=r"^Abc"),
    "OBP":  dict(pfam=["PF01395"], go=["GO:0005549"], name_regex=r"obp"),
    "P450": dict(pfam=["PF00067"], go=["GO:0020037", "GO:0016705"], name_regex=r"Cyp[0-9]"),
    "GST":  dict(pfam=["PF00043", "PF02798", "PF13409", "PF13410", "PF14497"],
                 go=["GO:0004364"], name_regex=r"^Gst"),
    "COE":  dict(pfam=["PF00135"], go=["GO:0004091"], name_regex=r"^Coe"),
}

_names = (pd.concat([gam[["GeneID", "GeneName"]], col[["GeneID", "GeneName"]]])
          .dropna(subset=["GeneID"]).drop_duplicates("GeneID"))
_ann = ann_by_gene.rename(columns={"gene_id": "GeneID"}).merge(_names, on="GeneID", how="outer")

FAM = {}
for fam, spec in FAMILIES.items():
    pfam_hit = _ann["pfam"].fillna("").apply(lambda s: any(p in s for p in spec["pfam"]))
    go_hit = _ann["go"].fillna("").apply(lambda s: any(g in s for g in spec["go"]))
    name_hit = _ann["GeneName"].fillna("").str.contains(spec["name_regex"], case=False, regex=True)
    FAM[fam] = set(_ann.loc[pfam_hit | go_hit | name_hit, "GeneID"])

# ── Gene coordinates ────────────────────────────────────────────────────
_GENES_CACHE = f"{SCRATCH}/genes.csv"
if not os.path.exists(_GENES_CACHE):
    import malariagen_data
    _gs = malariagen_data.Ag3(results_cache="~/ag3_cache", pre=True).geneset()
    (_gs[_gs["type"] == "gene"][["contig", "start", "end", "ID", "Name", "description"]]
     .rename(columns={"ID": "GeneID"}).to_csv(_GENES_CACHE, index=False))
genes_df = pd.read_csv(_GENES_CACHE)

# ── Sweep overlap (selection atlas H12, West African gamb/colu cohorts) ──
OFFSETS = {"2R": 0, "2L": 61545105, "3R": 0, "3L": 53200684, "X": 0}
ARM_TO_CONTIG = {"2R": "2RL", "2L": "2RL", "3R": "3RL", "3L": "3RL", "X": "X"}
WA = {"CI", "BF", "GH", "ML", "TG", "BJ", "GN", "GW", "GM", "SN", "SL", "LR", "NG", "NE", "CM"}

genes_df = genes_df.assign(
    atlas_contig=genes_df["contig"].map(ARM_TO_CONTIG),
    atlas_start=genes_df["start"] + genes_df["contig"].map(OFFSETS),
    atlas_end=genes_df["end"] + genes_df["contig"].map(OFFSETS),
).dropna(subset=["atlas_contig"])

signals = pd.read_csv("skills/selection-atlas/h12-signal-detection-all.csv")
signals["taxon"] = signals["cohort_id"].str.extract(r"_(gamb|colu|arab|bissau|fun)_")
signals["country"] = signals["cohort_id"].str[:2]
sig_wa = signals[signals["taxon"].isin(["gamb", "colu"]) & signals["country"].isin(WA)].copy()


def sweep_hits(row):
    """Signals whose span1 interval overlaps this gene."""
    s = sig_wa[sig_wa["contig"] == row["atlas_contig"]]
    return s[(s["span1_pstart"] <= row["atlas_end"]) & (s["span1_pstop"] >= row["atlas_start"])]


_OVERLAP_CACHE = f"{SCRATCH}/gene_sweeps.csv"
if os.path.exists(_OVERLAP_CACHE):
    _ov = pd.read_csv(_OVERLAP_CACHE)
else:
    recs = []
    for _, r in genes_df.iterrows():
        h = sweep_hits(r)
        recs.append(dict(
            GeneID=r["GeneID"], n_sweeps=len(h),
            max_delta_i=h["delta_i"].max() if len(h) else np.nan,
            max_stat_max=h["statistic_max"].max() if len(h) and "statistic_max" in h else np.nan,
            countries=",".join(sorted(h["country"].unique())) if len(h) else "",
            cohorts="; ".join(h["cohort_id"].unique()) if len(h) else "",
        ))
    _ov = pd.DataFrame(recs)
    _ov.to_csv(_OVERLAP_CACHE, index=False)

genes_df = genes_df.merge(_ov, on="GeneID", how="left")
genes_df["sweep_overlap"] = genes_df["n_sweeps"].fillna(0) > 0

# ── Convenience DE gene sets ────────────────────────────────────────────
def sig_up(de, ids):
    return set(de.loc[(de["padj"] < 0.05) & (de["log2FoldChange"] > 0) & de["GeneID"].isin(ids), "GeneID"])


abc_gam_de = sig_up(gam, FAM["ABC"])
abc_col_de = sig_up(col, FAM["ABC"])
obp_gam_de = sig_up(gam, FAM["OBP"])
obp_col_de = sig_up(col, FAM["OBP"])

if __name__ == "__main__":
    print({k: len(v) for k, v in FAM.items()})
    print("genes:", len(genes_df), "sweep-overlapping:", int(genes_df["sweep_overlap"].sum()))
    print("sweeps used:", len(sig_wa))
    print("DE sets:", len(abc_gam_de), len(abc_col_de), len(obp_gam_de), len(obp_col_de))
