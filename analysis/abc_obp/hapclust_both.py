"""ABCH1 haplotype clustering and NJ tree for BOTH species.

The selection-atlas signal at ABCH1 is present in 8 West African An. coluzzii cohorts and
0 of 16 An. gambiae cohorts, yet ABCH1 is more strongly upregulated in An. gambiae in Bouake.
Including both taxa in one tree shows directly whether the swept clade is coluzzii-specific.
"""
import os
import numpy as np
import allel
import anjl
from scipy.cluster.hierarchy import linkage, fcluster
import malariagen_data
from pathlib import Path

from bouake.ag3_sets import cap_per_cohort, unrestricted_sample_sets  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "abc_obp"

REGION = "2R:24,825,285-24,848,323"      # ABCH1
ANALYSIS = "gamb_colu"
OUT = str(RESULTS / "hapclust_both")
COUNTRIES = ["Burkina Faso", "Mali", "Ghana", "Cote d'Ivoire"]
MAX_PER_COHORT = 40
THRESH = 0.005

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
SAMPLE_SETS = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii", "gambiae"))

# Country names must match the metadata exactly: a near-miss spelling silently drops a whole
# country rather than raising, so check membership before querying.
available = set(ag3.sample_metadata(sample_sets=SAMPLE_SETS)["country"].unique())
missing = [c for c in COUNTRIES if c not in available]
assert not missing, f"country names absent from the metadata: {missing}"

base = "taxon in ['coluzzii', 'gambiae'] and country in {!r}".format(COUNTRIES)
meta_all = ag3.sample_metadata(sample_query=base, sample_sets=SAMPLE_SETS)
keep = cap_per_cohort(metadata=meta_all, cohort_column="cohort_admin2_year", max_per_cohort=MAX_PER_COHORT)
print("samples kept:", len(keep), flush=True)
print(keep.groupby(["country", "taxon"]).size().to_string(), flush=True)

query = "sample_id in {!r}".format(list(keep["sample_id"]))
ds = ag3.haplotypes(region=REGION, analysis=ANALYSIS, sample_query=query, sample_sets=SAMPLE_SETS)
gt = allel.GenotypeDaskArray(ds["call_genotype"].data)
ht = np.asarray(gt.to_haplotypes().compute())
samples = ds["sample_id"].values
X = ht.T.astype(np.float32)
n_hap, n_sites = X.shape
print(f"haplotypes {n_hap} x sites {n_sites}", flush=True)

meta = ag3.sample_metadata(sample_query=query, sample_sets=SAMPLE_SETS).set_index("sample_id").loc[samples].reset_index()
hap_meta = meta.loc[np.repeat(meta.index, 2)].reset_index(drop=True)

rs = X.sum(axis=1)
D = (rs[:, None] + rs[None, :] - 2 * (X @ X.T)) / n_sites
np.fill_diagonal(D, 0.0)
D = np.clip(D, 0, None).astype(np.float32)

Z = linkage(D[np.triu_indices(n_hap, 1)], method="complete")
hap_meta["cluster"] = fcluster(Z, t=THRESH, criterion="distance")
sizes = hap_meta["cluster"].value_counts()
major = sizes[sizes >= 20].index.tolist()
print(f"\n{len(major)} clusters with n>=20; largest {sizes.iloc[0]}", flush=True)
print(hap_meta[hap_meta["cluster"].isin(major)]
      .groupby(["cluster", "taxon"]).size().unstack(fill_value=0)
      .sort_values(by=list(hap_meta.taxon.unique())[0], ascending=False).head(12).to_string(),
      flush=True)

tree = anjl.canonical_nj(D)
np.save(f"{OUT}/njt_Z.npy", np.asarray(tree))
np.save(f"{OUT}/linkage.npy", Z)
hap_meta.to_csv(f"{OUT}/hap_meta.csv", index=False)
print("\ndone", flush=True)
