"""A5. The top-N ABC sweep enrichment is inflated by tandem arrays (4 adjacent ABCCs share one
signal). Re-test collapsing genes within 200 kb on the same arm into a single independent locus.
Also identify the non-canonical 'ABC' genes pulled in by the broad PFAM/GO definition."""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from setup import gam, genes_df, abc_gam_de, N_PERM

pd.set_option("display.width", 220)
CLUSTER_BP = 200_000


def collapse(df):
    """Group genes on the same arm within CLUSTER_BP into loci; a locus is 'swept' if any member is."""
    out = []
    for arm, sub in df.sort_values(["contig", "start"]).groupby("contig"):
        locus, last = 0, None
        for _, r in sub.iterrows():
            if last is not None and r["start"] - last > CLUSTER_BP:
                locus += 1
            last = r["start"]
            out.append((r["GeneID"], f"{arm}:{locus}"))
    m = dict(out)
    d = df.assign(locus=df["GeneID"].map(m))
    return (d.groupby("locus")
            .agg(contig=("contig", "first"), n_genes=("GeneID", "size"),
                 swept=("sweep_overlap", "any"),
                 genes=("GeneID", lambda x: ",".join(x)))
            .reset_index())


# Background: collapse the whole genome the same way
# NOTE: proximity-chaining the whole genome merges everything (genes are denser than CLUSTER_BP),
# so the background stays at gene resolution. Only the OBSERVED set is collapsed: random draws of
# scattered genes are almost never tandem, so the null needs no collapsing.
print(f"Background: {len(genes_df)} genes, sweep rate {genes_df['sweep_overlap'].mean():.3f}")

arm_pool = {a: genes_df.loc[genes_df["contig"] == a, "sweep_overlap"].to_numpy()
            for a in genes_df["contig"].unique()}

abc = (gam[gam["GeneID"].isin(abc_gam_de)][["GeneID", "GeneName", "log2FoldChange"]]
       .merge(genes_df[["GeneID", "contig", "start", "sweep_overlap", "n_sweeps", "Name"]],
              on="GeneID", how="left")
       .sort_values("log2FoldChange", ascending=False).reset_index(drop=True))
abc["label"] = abc["GeneName"].fillna(abc["Name"]).fillna(abc["GeneID"])


def perm_loci(sub, n_perm=N_PERM, seed=42):
    loci = collapse(sub[["GeneID", "contig", "start", "sweep_overlap"]])
    obs_hits = int(loci["swept"].sum())
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm, dtype=np.int32)
    for arm, n in loci["contig"].value_counts().items():
        pool = arm_pool[arm]
        idx = np.argsort(rng.random((n_perm, len(pool))), axis=1)[:, :n]
        null += pool[idx].sum(axis=1)
    bg_hits = int(genes_df["sweep_overlap"].sum()) - obs_hits
    bg_n = len(genes_df) - len(loci)
    orr, fp = fisher_exact([[obs_hits, len(loci) - obs_hits], [bg_hits, bg_n - bg_hits]],
                           alternative="greater")
    return dict(n_loci=len(loci), hits=obs_hits, rate=obs_hits / len(loci),
                null_mean=null.mean(), perm_p=float((null >= obs_hits).mean()),
                OR=orr, fisher_p=fp)


print("\nA5. Top-N ABC sweep enrichment AFTER collapsing tandem arrays into independent loci")
rows = []
for n in [5, 10, 15, 20, 25, 30]:
    rows.append(dict(top_n=n, **perm_loci(abc.head(n))))
print(pd.DataFrame(rows).round(4).to_string(index=False))

print("\n  Loci making up the top 20 ABCs:")
print(collapse(abc.head(20)[["GeneID", "contig", "start", "sweep_overlap"]])
      .assign(labels=lambda d: d["genes"].apply(
          lambda g: ",".join(abc.set_index("GeneID").loc[g.split(","), "label"])))
      [["locus", "n_genes", "swept", "labels"]].to_string(index=False))

# ── Which 'ABCs' are not canonical ABC transporters? ─────────────────────
print("\nA6. Genes in the ABC set whose name is not Abc* — check the family definition is not "
      "pulling in unrelated ATP-binding proteins")
odd = abc[~abc["label"].str.contains("ABC", case=False, na=False)]
print(odd[["label", "GeneID", "contig", "log2FoldChange", "n_sweeps"]].to_string(index=False))
