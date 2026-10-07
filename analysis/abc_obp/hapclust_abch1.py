"""Haplotype clustering at ABCH1.

Answers two questions the H12 traces cannot:
  1. Is the signal one sweep or several? -> how many large, low-diversity haplotype clusters.
  2. What SNPs tag each swept cluster? -> candidates for tracking frequency through time.

Scoped to the ABCH1 gene itself. A wider window pulls in recombinants from the flanks and
blurs the cluster structure, as well as costing more sites.

Hamming distances are computed by matrix multiplication rather than allel.pairwise_distance:
for binary haplotypes, mismatches = a.sum + b.sum - 2*(a·b), so the whole matrix is one dot
product. That is the difference between seconds and tens of minutes at this sample size.

Run on hyperion. Writes to ./hapclust_out/.
"""

import os
import numpy as np
import pandas as pd
import allel
from scipy.cluster.hierarchy import linkage, fcluster
import malariagen_data
from pathlib import Path

from bouake.ag3_sets import cap_per_cohort, unrestricted_sample_sets  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "abc_obp"

REGION = "2R:24,825,285-24,848,323"    # ABCH1 only
ANALYSIS = "gamb_colu"
OUT = str(RESULTS / "hapclust")
COUNTRIES = ["Burkina Faso", "Mali", "Ghana", "Cote d'Ivoire"]
MAX_PER_COHORT = 50                    # keep cohorts comparable and the matrix tractable
THRESH = 0.005                         # distance cut defining a near-identical haplotype set

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
SAMPLE_SETS = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii",))

# ── Stratified subsample so no single cohort dominates ──────────────────
base_query = "taxon == 'coluzzii' and country in {!r}".format(COUNTRIES)
meta_all = ag3.sample_metadata(sample_query=base_query, sample_sets=SAMPLE_SETS)
keep = cap_per_cohort(metadata=meta_all, cohort_column="cohort_admin2_year", max_per_cohort=MAX_PER_COHORT)
print(f"samples: {len(meta_all)} -> {len(keep)} after capping at {MAX_PER_COHORT}/cohort",
      flush=True)

query = "sample_id in {!r}".format(list(keep["sample_id"]))
ds = ag3.haplotypes(region=REGION, analysis=ANALYSIS, sample_query=query, sample_sets=SAMPLE_SETS)
print("dataset:", dict(ds.sizes), flush=True)

gt = allel.GenotypeDaskArray(ds["call_genotype"].data)
ht = np.asarray(gt.to_haplotypes().compute())      # sites x haplotypes
pos = ds["variant_position"].values
samples = ds["sample_id"].values
X = ht.T.astype(np.float32)                        # haplotypes x sites
n_hap, n_sites = X.shape
print(f"haplotypes: {n_hap} over {n_sites} sites", flush=True)

meta = ag3.sample_metadata(sample_query=query, sample_sets=SAMPLE_SETS).set_index("sample_id").loc[samples].reset_index()
hap_meta = meta.loc[np.repeat(meta.index, 2)].reset_index(drop=True)
hap_meta["hap_id"] = [f"{s}_{i % 2}" for i, s in enumerate(np.repeat(samples, 2))]

# ── Pairwise Hamming via one dot product ────────────────────────────────
rs = X.sum(axis=1)
inner = X @ X.T
D = (rs[:, None] + rs[None, :] - 2 * inner) / n_sites
np.fill_diagonal(D, 0.0)
D = np.clip(D, 0, None)
condensed = D[np.triu_indices(n_hap, 1)]
print(f"distance matrix done; median pairwise {np.median(condensed):.4f}", flush=True)

Z = linkage(condensed, method="complete")
for t in (0.001, 0.002, 0.005, 0.01, 0.02):
    cl = fcluster(Z, t=t, criterion="distance")
    sizes = pd.Series(cl).value_counts()
    big = sizes[sizes >= 20]
    print(f"threshold {t}: {len(sizes)} clusters, {len(big)} with n>=20 "
          f"covering {big.sum()}/{len(cl)} ({100*big.sum()/len(cl):.1f}%)", flush=True)

hap_meta["cluster"] = fcluster(Z, t=THRESH, criterion="distance")
sizes = hap_meta["cluster"].value_counts()
major = sizes[sizes >= 20].index.tolist()
print(f"\n{len(major)} major clusters (n>=20) at threshold {THRESH}", flush=True)
print(hap_meta[hap_meta["cluster"].isin(major)]
      .groupby(["cluster", "country"]).size().unstack(fill_value=0).to_string(), flush=True)
print("\nBurkina Faso by year:", flush=True)
print(hap_meta[hap_meta["cluster"].isin(major) & (hap_meta["country"] == "Burkina Faso")]
      .groupby(["cluster", "year"]).size().unstack(fill_value=0).to_string(), flush=True)

print("\nwithin-cluster mean pairwise distance (hard sweep -> ~0):", flush=True)
for c in major:
    idx = np.where(hap_meta["cluster"].values == c)[0]
    sub = D[np.ix_(idx, idx)]
    print(f"  cluster {c}: n={len(idx)}, mean={sub[np.triu_indices(len(idx),1)].mean():.5f}",
          flush=True)

# ── Tagging SNPs: near-fixed inside a cluster, rare outside ─────────────
ref = ds["variant_allele"].values[:, 0].astype(str)
alt = ds["variant_allele"].values[:, 1].astype(str)
rows = []
for c in major:
    inside = hap_meta["cluster"].values == c
    f_in, f_out = X[inside].mean(axis=0), X[~inside].mean(axis=0)
    delta = f_in - f_out
    for i in np.argsort(-delta)[:40]:
        rows.append(dict(cluster=c, n_hap=int(inside.sum()), position=int(pos[i]),
                         ref=ref[i], alt=alt[i], freq_in=round(float(f_in[i]), 4),
                         freq_out=round(float(f_out[i]), 4), delta=round(float(delta[i]), 4)))
tags = pd.DataFrame(rows)
tags.to_csv(f"{OUT}/tagging_snps.csv", index=False)
print("\ntop tagging SNPs per cluster (delta > 0.5):", flush=True)
print(tags[tags["delta"] > 0.5].groupby("cluster").head(8).to_string(index=False), flush=True)

hap_meta.to_csv(f"{OUT}/hap_clusters.csv", index=False)
np.save(f"{OUT}/linkage.npy", Z)
np.save(f"{OUT}/positions.npy", pos)
print("\ndone", flush=True)
