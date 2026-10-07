"""B. T3-redux (clustering), T4 (cross-study rank shift), and a background-matched
species-divergence test."""

import os
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from setup import (gam, col, FAM, genes_df, N_PERM, SCRATCH,
                   abc_gam_de, abc_col_de, obp_gam_de, obp_col_de)

pd.set_option("display.width", 240)
rng_global = np.random.default_rng(42)

# =======================================================================
# T3-redux — clustering with nearest-neighbour + arm-concentration
# =======================================================================
GENOME_PENALTY = 3e8


def nn_dist(df):
    """Median nearest-neighbour distance (bp); genes with no same-arm sibling take a penalty."""
    if len(df) < 2:
        return np.nan
    pos = df[["contig", "start"]].to_numpy()
    nn = []
    for i in range(len(pos)):
        same = [abs(pos[j, 1] - pos[i, 1]) for j in range(len(pos))
                if j != i and pos[j, 0] == pos[i, 0]]
        nn.append(min(same) if same else GENOME_PENALTY)
    return float(np.median(nn))


def arm_conc(df):
    return float(df["contig"].value_counts(normalize=True).iloc[0]) if len(df) else np.nan


def cluster_test(family_ids, de_ids, label, n_perm=N_PERM, seed=42):
    fam = genes_df[genes_df["GeneID"].isin(family_ids)][["GeneID", "contig", "start"]].drop_duplicates("GeneID")
    de = fam[fam["GeneID"].isin(de_ids)]
    if len(de) < 3:
        return None
    obs_nn, obs_arm = nn_dist(de), arm_conc(de)
    rng = np.random.default_rng(seed)
    pool = fam["GeneID"].to_numpy()
    null_nn = np.empty(n_perm)
    null_arm = np.empty(n_perm)
    fam_idx = fam.set_index("GeneID")
    for i in range(n_perm):
        pick = rng.choice(pool, size=len(de), replace=False)
        d = fam_idx.loc[pick].reset_index()
        null_nn[i] = nn_dist(d)
        null_arm[i] = arm_conc(d)
    return dict(set=label, n_family=len(fam), n_de=len(de),
                obs_nn_kb=obs_nn / 1e3, null_nn_kb=np.median(null_nn) / 1e3,
                p_nn=float(np.mean(null_nn <= obs_nn)),
                obs_arm_conc=obs_arm, null_arm_conc=np.median(null_arm),
                p_arm=float(np.mean(null_arm >= obs_arm)),
                arms=str(de["contig"].value_counts().to_dict()))


print("=" * 110)
print("T3-redux. Genomic clustering of the DE response (nearest-neighbour + arm concentration)")
print("=" * 110)
rows = [r for r in [
    cluster_test(FAM["OBP"], obp_col_de, "coluzzii OBP-up"),
    cluster_test(FAM["OBP"], obp_gam_de, "gambiae OBP-up"),
    cluster_test(FAM["ABC"], abc_gam_de, "gambiae ABC-up"),
    cluster_test(FAM["ABC"], abc_col_de, "coluzzii ABC-up"),
] if r]
t3 = pd.DataFrame(rows)
for c in ["p_nn", "p_arm"]:
    t3[c + "_adj"] = multipletests(t3[c], method="fdr_bh")[1]
print(t3.round(4).to_string(index=False))

# =======================================================================
# T4 — cross-study rank shift (AnoExpress meta fold change)
# =======================================================================
import anoexpress as ae

CACHE = f"{SCRATCH}/ae_meta.csv"
if os.path.exists(CACHE):
    meta = pd.read_csv(CACHE)
else:
    cand = ae.load_candidates(analysis="gamb_colu", func=np.nanmedian).reset_index()
    cand = cand.loc[:, ~cand.columns.duplicated()]
    fc = [c for c in cand.columns if "log2" in c.lower()][0]
    meta = cand[["GeneID", fc]].rename(columns={fc: "ae_median_log2FC"}).drop_duplicates("GeneID")
    meta.to_csv(CACHE, index=False)
print("\nAnoExpress meta columns:", list(meta.columns)[:8], "| n =", len(meta))

FC_COL = [c for c in meta.columns if "log2" in c.lower()][0]


def t4(de, fam_ids, label):
    d = de[de["GeneID"].isin(fam_ids)].merge(meta[["GeneID", FC_COL]], on="GeneID", how="inner").dropna(subset=[FC_COL])
    up = d[(d["padj"] < 0.05) & (d["log2FoldChange"] > 0)][FC_COL]
    ns = d[~((d["padj"] < 0.05) & (d["log2FoldChange"] > 0))][FC_COL]
    if len(up) < 3 or len(ns) < 3:
        return None
    u, p = stats.mannwhitneyu(up, ns, alternative="greater")
    # rank-biserial effect size
    rb = 2 * u / (len(up) * len(ns)) - 1
    return dict(set=label, n_up=len(up), n_ns=len(ns),
                med_meta_up=up.median(), med_meta_ns=ns.median(),
                rank_biserial=rb, U=u, p=p)


print("\n" + "=" * 110)
print("T4. Cross-study rank shift — is Bouaké-sig-up associated with a higher AnoExpress meta log2FC?")
print("   (one-sided: Bouaké-up > Bouaké-ns; a null result means the Bouaké response is NOT recovered")
print("    by the pyrethroid-dominated meta-analysis, i.e. it is genuinely site/insecticide-specific)")
print("=" * 110)
rows = [r for r in [
    t4(gam, FAM["ABC"], "ABC gambiae"), t4(col, FAM["ABC"], "ABC coluzzii"),
    t4(gam, FAM["OBP"], "OBP gambiae"), t4(col, FAM["OBP"], "OBP coluzzii"),
    t4(gam, FAM["P450"], "P450 gambiae"), t4(col, FAM["P450"], "P450 coluzzii"),
] if r]
t4df = pd.DataFrame(rows)
t4df["padj"] = multipletests(t4df["p"], method="fdr_bh")[1]
print(t4df.round(4).to_string(index=False))

# Threshold comparison the user asked for: shared/novel counts under each cutoff
print("\nT4b. 'Shared candidate' counts under each candidate threshold")
for label, de, fam_ids in [("ABC gambiae", gam, FAM["ABC"]), ("OBP coluzzii", col, FAM["OBP"])]:
    d = de[de["GeneID"].isin(fam_ids)].merge(meta[["GeneID", FC_COL]], on="GeneID", how="left")
    up = d[(d["padj"] < 0.05) & (d["log2FoldChange"] > 0)]
    n_in = up[FC_COL].notna().sum()
    print(f"  {label}: {len(up)} sig-up, {n_in} present in AnoExpress")
    for thr in [0.3, 0.5]:
        print(f"     meta log2FC > {thr}: shared={int((up[FC_COL] > thr).sum())}, "
              f"novel={int((up[FC_COL] <= thr).sum())}")

# =======================================================================
# Species divergence, background-matched
# =======================================================================
paired = (gam[["GeneID", "log2FoldChange", "padj"]].rename(columns={"log2FoldChange": "l2fc_gam", "padj": "padj_gam"})
          .merge(col[["GeneID", "log2FoldChange", "padj"]].rename(columns={"log2FoldChange": "l2fc_col", "padj": "padj_col"}),
                 on="GeneID", how="inner")
          .dropna(subset=["l2fc_gam", "l2fc_col"]))
paired["delta"] = paired["l2fc_gam"] - paired["l2fc_col"]
paired = paired.merge(genes_df[["GeneID", "contig"]], on="GeneID", how="inner")
print(f"\nBackground: {len(paired)} genes tested in both species; "
      f"genome-wide mean delta (gam - col) = {paired['delta'].mean():.4f}")

arm_delta = {a: paired.loc[paired["contig"] == a, "delta"].to_numpy() for a in paired["contig"].unique()}


def divergence_test(fam_ids, label, n_perm=N_PERM, seed=42):
    """Is the family's gambiae-minus-coluzzii response shift larger than an arm-matched random set?
    Two-sided: controls for the genome-wide species offset by construction (null uses same statistic)."""
    obs_df = paired[paired["GeneID"].isin(fam_ids)]
    if len(obs_df) < 10:
        return None
    obs = obs_df["delta"].mean()
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm)
    total = 0
    for arm, n in obs_df["contig"].value_counts().items():
        pool = arm_delta[arm]
        idx = np.argsort(rng.random((n_perm, len(pool))), axis=1)[:, :n]
        null += pool[idx].sum(axis=1)
        total += n
    null /= total
    p_two = float(np.mean(np.abs(null - null.mean()) >= abs(obs - null.mean())))
    return dict(set=label, n=len(obs_df), obs_mean_delta=obs,
                null_mean=null.mean(), null_sd=null.std(),
                z=(obs - null.mean()) / null.std(), p=p_two)


print("\n" + "=" * 110)
print("Species divergence, arm-matched permutation — mean(log2FC_gambiae - log2FC_coluzzii) per family")
print("  vs random gene sets matched on chromosome-arm distribution. This is the test T1's Wilcoxon")
print("  lacked: it controls for the genome-wide difference in response magnitude between species.")
print("=" * 110)
rows = [r for r in [divergence_test(FAM[f], f) for f in ["ABC", "OBP", "P450", "GST", "COE"]] if r]
dv = pd.DataFrame(rows)
dv["padj"] = multipletests(dv["p"], method="fdr_bh")[1]
print(dv.round(4).to_string(index=False))

t3.to_csv("results/abc_obp/t3_redux.csv", index=False)
t4df.to_csv("results/abc_obp/t4.csv", index=False)
dv.to_csv("results/abc_obp/divergence.csv", index=False)
