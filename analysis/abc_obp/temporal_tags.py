"""Temporal frequency of the cluster-168 tagging SNPs in Burkina Faso.

The cluster frequencies come from a subsample capped at 50 samples per cohort. This recomputes
allele frequencies at the top tagging SNPs using every available sample, by province and year,
so the trajectory is not limited by that cap.

Run on hyperion. Writes to ./temporal_out/.
"""

import os
import numpy as np
import pandas as pd
import allel
import malariagen_data
from pathlib import Path

from bouake.ag3_sets import unrestricted_sample_sets  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "abc_obp"

REGION = "2R:24,825,000-24,848,500"
OUT = str(RESULTS / "temporal")
# Top tags for cluster 168 (freq_in ~1.0, freq_out < 0.17)
TAGS = [24840512, 24840627, 24840626, 24831203, 24830530, 24836045, 24836042,
        24841324, 24839975, 24827017, 24827562, 24846834, 24835989, 24837992, 24844319]

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
SAMPLE_SETS = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii",))

query = ("taxon == 'coluzzii' and country in "
         "['Burkina Faso', 'Mali', 'Ghana', \"Cote d'Ivoire\"]")
ds = ag3.snp_calls(region=REGION, sample_query=query, sample_sets=SAMPLE_SETS, site_mask="gamb_colu")
pos = ds["variant_position"].values
keep = np.isin(pos, TAGS)
print(f"tags found in call set: {keep.sum()} / {len(TAGS)}", flush=True)

gt = allel.GenotypeDaskArray(ds["call_genotype"].data[keep]).compute()
meta = ag3.sample_metadata(sample_query=query, sample_sets=SAMPLE_SETS).set_index("sample_id") \
    .loc[ds["sample_id"].values].reset_index()
tag_pos = pos[keep]

rows = []
for (country, admin2, year), grp in meta.groupby(["country", "admin2_name", "year"]):
    if len(grp) < 10:
        continue
    idx = grp.index.values
    sub = gt.take(idx, axis=1)
    ac = sub.count_alleles()
    for j, p in enumerate(tag_pos):
        n = ac[j].sum()
        if n == 0:
            continue
        rows.append(dict(country=country, admin2=admin2, year=int(year), position=int(p),
                         n_samples=len(grp), an=int(n),
                         alt_freq=float(ac[j, 1:].sum() / n)))
df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/tag_freqs_by_cohort.csv", index=False)

print("\nBurkina Faso, mean tag frequency by province and year:", flush=True)
bf = df[df["country"] == "Burkina Faso"]
piv = bf.groupby(["admin2", "year"]).agg(mean_freq=("alt_freq", "mean"),
                                         n_samples=("n_samples", "first")).reset_index()
print(piv.to_string(index=False), flush=True)

print("\nPer-tag trajectory, Houet:", flush=True)
ho = bf[bf["admin2"] == "Houet"].pivot_table(index="year", columns="position",
                                             values="alt_freq")
print(ho.round(3).to_string(), flush=True)
print("\ndone", flush=True)
