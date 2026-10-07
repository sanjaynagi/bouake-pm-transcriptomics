"""A. Sweep overlap focused on the MOST upregulated ABCs, rather than all sig-up ABCs."""

import numpy as np
import pandas as pd
from scipy import stats
from setup import (gam, col, FAM, genes_df, sig_wa, abc_gam_de, obp_col_de, N_PERM)

pd.set_option("display.width", 250)
pd.set_option("display.max_colwidth", 70)

coords = genes_df.set_index("GeneID")

# ── A1. Per-gene detail, ABC sig-up in gambiae, ranked by log2FC ─────────
abc = (gam[gam["GeneID"].isin(abc_gam_de)][["GeneID", "GeneName", "log2FoldChange", "padj"]]
       .merge(genes_df[["GeneID", "contig", "start", "end", "Name",
                        "n_sweeps", "max_delta_i", "countries", "cohorts"]],
              on="GeneID", how="left")
       .sort_values("log2FoldChange", ascending=False)
       .reset_index(drop=True))
abc["label"] = abc["GeneName"].fillna(abc["Name"]).fillna(abc["GeneID"])
abc["rank"] = np.arange(1, len(abc) + 1)

print("=" * 100)
print("A1. Bouaké gambiae sig-up ABCs, ranked by log2FC, with West-African H12 sweep overlap")
print("=" * 100)
print(abc[["rank", "label", "GeneID", "contig", "log2FoldChange", "padj",
           "n_sweeps", "max_delta_i", "countries"]].to_string(index=False))

# ── A2. Do sweep-overlapping ABCs have higher log2FC? ────────────────────
sw = abc[abc["n_sweeps"] > 0]["log2FoldChange"]
no = abc[abc["n_sweeps"] == 0]["log2FoldChange"]
u, p = stats.mannwhitneyu(sw, no, alternative="two-sided")
print(f"\nA2. log2FC of swept vs non-swept sig-up ABCs: "
      f"median {sw.median():.3f} (n={len(sw)}) vs {no.median():.3f} (n={len(no)}); "
      f"Mann-Whitney U={u:.0f}, p={p:.4f}")

# Rank-correlation across ALL ABCs (not just sig-up): is n_sweeps related to log2FC?
abc_all = (gam[gam["GeneID"].isin(FAM["ABC"])][["GeneID", "GeneName", "log2FoldChange", "padj"]]
           .merge(genes_df[["GeneID", "contig", "n_sweeps"]], on="GeneID", how="inner")
           .dropna(subset=["log2FoldChange", "n_sweeps"]))
rho, prho = stats.spearmanr(abc_all["log2FoldChange"], abc_all["n_sweeps"])
print(f"    Spearman(log2FC, n_sweeps) across all {len(abc_all)} ABCs: rho={rho:.3f}, p={prho:.4f}")

# ── A3. Arm-matched permutation, restricted to the TOP-N ABCs by log2FC ──
arm_pool = {a: genes_df.loc[genes_df["contig"] == a, "sweep_overlap"].to_numpy()
            for a in genes_df["contig"].unique()}


def perm_p(gene_ids, n_perm=N_PERM, seed=42):
    obs = genes_df[genes_df["GeneID"].isin(set(gene_ids))]
    if obs.empty:
        return None
    obs_hits = int(obs["sweep_overlap"].sum())
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm, dtype=np.int32)
    for arm, n in obs["contig"].value_counts().items():
        pool = arm_pool[arm]
        idx = np.argsort(rng.random((n_perm, len(pool))), axis=1)[:, :n]
        null += pool[idx].sum(axis=1)
    from scipy.stats import fisher_exact
    bg_hits = int(genes_df["sweep_overlap"].sum()) - obs_hits
    bg_n = len(genes_df) - len(obs)
    orr, fp = fisher_exact([[obs_hits, len(obs) - obs_hits], [bg_hits, bg_n - bg_hits]],
                           alternative="greater")
    return dict(n=len(obs), hits=obs_hits, rate=obs_hits / len(obs),
                null_mean=null.mean(), null_sd=null.std(),
                perm_p=float((null >= obs_hits).mean()), OR=orr, fisher_p=fp)


print("\nA3. Arm-matched permutation enrichment as a function of how many top ABCs are included")
rows = []
for n in [5, 10, 15, 20, 25, 30]:
    r = perm_p(abc.head(n)["GeneID"])
    rows.append(dict(top_n=n, min_log2FC=abc.head(n)["log2FoldChange"].min(), **r))
print(pd.DataFrame(rows).round(4).to_string(index=False))

# Same for OBPs (coluzzii) as comparator
obp = (col[col["GeneID"].isin(obp_col_de)][["GeneID", "GeneName", "log2FoldChange", "padj"]]
       .merge(genes_df[["GeneID", "contig", "n_sweeps"]], on="GeneID", how="left")
       .sort_values("log2FoldChange", ascending=False))
print("\n    OBP (coluzzii) comparator:")
rows = []
for n in [5, 10, 16]:
    r = perm_p(obp.head(n)["GeneID"])
    rows.append(dict(top_n=n, **r))
print(pd.DataFrame(rows).round(4).to_string(index=False))

# ── A4. Detail on the swept ABCs — which cohorts, how strong ─────────────
print("\n" + "=" * 100)
print("A4. Sweep detail for every sig-up ABC that overlaps a signal")
print("=" * 100)
for _, r in abc[abc["n_sweeps"] > 0].iterrows():
    g = coords.loc[r["GeneID"]]
    h = sig_wa[(sig_wa["contig"] == g["atlas_contig"]) &
               (sig_wa["span1_pstart"] <= g["atlas_end"]) &
               (sig_wa["span1_pstop"] >= g["atlas_start"])]
    print(f"\n--- {r['label']} ({r['GeneID']}) {r['contig']}:{int(g['start']):,}-{int(g['end']):,} "
          f"| log2FC {r['log2FoldChange']:.2f}, padj {r['padj']:.2e} | {len(h)} signals")
    cols = [c for c in ["cohort_id", "country", "taxon", "span1_pstart", "span1_pstop",
                        "delta_i", "statistic_max", "focus_pstart", "focus_pstop"] if c in h.columns]
    print(h[cols].sort_values("delta_i", ascending=False).head(12).to_string(index=False))

abc.to_csv("results/abc_obp/abc_gam_sweeps.csv", index=False)
