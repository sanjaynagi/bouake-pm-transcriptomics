"""Population branch statistic across 2R, An. coluzzii cohorts with An. arabiensis as outgroup.

PBS is not in malariagen_data (PR #1299 is unmerged), so it is computed here from Hudson Fst:

    T(X,Y) = -log(1 - Fst(X,Y))
    PBS(A) = (T(A,B) + T(A,C) - T(B,C)) / 2

with A the focal An. coluzzii cohort, B An. gambiae from the same country, and C An. arabiensis
from Burkina Faso. A branch specific to A is what distinguishes a coluzzii sweep from ancestral
differentiation shared with gambiae, which H12 alone cannot separate.

Sites use the gamb_colu_arab mask so all three taxa are called on the same accessible sites.
Fst is a ratio of averages within each window, so numerator and denominator are summed
separately rather than averaging per-site Fst.

Restricted to a 2.5 Mb window around ABCH1: percentiles below are therefore local
background, not genome-wide.

Run on hyperion. Writes to ./pbs_out/.
"""

import os
import numpy as np
import pandas as pd
import allel
import malariagen_data
from pathlib import Path

from bouake.ag3_sets import unrestricted_sample_sets  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "abc_obp"

CONTIG = "2R"
REGION_START, REGION_STOP = 23_500_000, 26_000_000   # the locus plus flanks for local context
SITE_MASK = "gamb_colu_arab"
OUT = str(RESULTS / "pbs")
BLOCK = 5_000_000          # bp loaded at a time, to bound memory
WINDOW = 10_000            # bp per PBS window; PBS is choppier than H12 at equal span
MIN_SNPS = 40              # segregating sites needed before a window is reported
CAP_COHORT = 50
CAP_OUTGROUP = 100
SEED = 42

COHORTS = [
    ("BF-02_Comoe_colu_2012",              "year",    "Burkina Faso"),
    ("BF-09_Houet_colu_2012_Q3",           "quarter", "Burkina Faso"),
    ("BF-02_Comoe_colu_2015",              "year",    "Burkina Faso"),
    ("BF-02_Comoe_colu_2016",              "year",    "Burkina Faso"),
    ("BF-09_Houet_colu_2014_Q3",           "quarter", "Burkina Faso"),
    ("BF-02_Comoe_colu_2011",              "year",    "Burkina Faso"),
    ("ML-3_Yanfolila_colu_2012_Q4",        "quarter", "Mali"),
    ("GH-AH_Adansi-Akrofuom_colu_2018_Q4", "quarter", "Ghana"),
]

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
SAMPLE_SETS = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii", "gambiae", "arabiensis"))
meta = ag3.sample_metadata(sample_sets=SAMPLE_SETS)
rng = np.random.RandomState(SEED)


def cap(df, n):
    return df.sample(min(len(df), n), random_state=SEED)


# ── Assemble the sample sets: focal cohorts, their gambiae comparators, the outgroup ──
groups = {}
for cohort_id, level, country in COHORTS:
    col = "cohort_admin2_quarter" if level == "quarter" else "cohort_admin2_year"
    groups[cohort_id] = cap(meta[meta[col] == cohort_id], CAP_COHORT)

for country in sorted({c for _, _, c in COHORTS}):
    sel = meta[(meta.taxon == "gambiae") & (meta.country == country)]
    groups[f"gamb_{country}"] = cap(sel, CAP_COHORT)

groups["arab_BF"] = cap(meta[(meta.taxon == "arabiensis")
                             & (meta.country == "Burkina Faso")], CAP_OUTGROUP)

for name, df in groups.items():
    print(f"{name:38s} n={len(df)}", flush=True)

all_ids = sorted({s for df in groups.values() for s in df.sample_id})
print(f"\ntotal samples: {len(all_ids)}", flush=True)

# ── One pass over the contig, accumulating per-window Fst numerator and denominator ──
span = REGION_STOP - REGION_START
pairs = []
for cohort_id, _, country in COHORTS:
    pairs += [(cohort_id, f"gamb_{country}"), (cohort_id, "arab_BF")]
pairs += [(f"gamb_{c}", "arab_BF") for c in sorted({c for _, _, c in COHORTS})]
pairs = sorted(set(pairs))

n_win = int(np.ceil(span / WINDOW))
num = {p: np.zeros(n_win) for p in pairs}
den = {p: np.zeros(n_win) for p in pairs}
cnt = np.zeros(n_win)

for start in range(REGION_START, REGION_STOP, BLOCK):
    stop = min(start + BLOCK - 1, REGION_STOP)
    region = f"{CONTIG}:{start:,}-{stop:,}"
    ds = ag3.snp_calls(region=region, site_mask=SITE_MASK,
                       sample_query="sample_id in {!r}".format(all_ids), sample_sets=SAMPLE_SETS)
    pos = ds["variant_position"].values
    ids = list(ds["sample_id"].values)
    index = {s: i for i, s in enumerate(ids)}
    subpops = {name: [index[s] for s in df.sample_id if s in index]
               for name, df in groups.items()}

    gt = allel.GenotypeDaskArray(ds["call_genotype"].data)
    acs = {k: np.asarray(v.compute())
           for k, v in gt.count_alleles_subpops(subpops, max_allele=3).items()}

    # Invariant sites contribute 0/0 to Hudson Fst, so drop them: same answer, less arithmetic
    total = sum(acs.values())
    segregating = (total > 0).sum(axis=1) > 1
    acs = {k: v[segregating] for k, v in acs.items()}
    pos = pos[segregating]

    win = np.clip((pos - REGION_START) // WINDOW, 0, n_win - 1)
    cnt += np.bincount(win, minlength=n_win)

    for a, b in pairs:
        n_, d_ = allel.hudson_fst(acs[a], acs[b])
        n_ = np.nan_to_num(np.asarray(n_))
        d_ = np.nan_to_num(np.asarray(d_))
        num[(a, b)] += np.bincount(win, weights=n_, minlength=n_win)
        den[(a, b)] += np.bincount(win, weights=d_, minlength=n_win)

    print(f"[{region}] {segregating.sum():,} segregating of {len(segregating):,} "
          f"accessible sites", flush=True)

# ── Fst -> divergence -> PBS ──────────────────────────────────────────────
def divergence(pair):
    with np.errstate(divide="ignore", invalid="ignore"):
        fst = np.where(den[pair] > 0, num[pair] / den[pair], np.nan)
    fst = np.clip(fst, 0, 0.9999)
    return -np.log(1 - fst)


mid = REGION_START + (np.arange(n_win) + 0.5) * WINDOW
enough = cnt >= MIN_SNPS
rows = {"position": mid, "n_segregating": cnt}
for cohort_id, _, country in COHORTS:
    b, c = f"gamb_{country}", "arab_BF"
    t_ab, t_ac, t_bc = divergence((cohort_id, b)), divergence((cohort_id, c)), divergence((b, c))
    pbs = (t_ab + t_ac - t_bc) / 2
    rows[cohort_id] = np.where(enough, pbs, np.nan)

df = pd.DataFrame(rows)
df.to_csv(f"{OUT}/pbs_2R.csv", index=False)

# ── Where does ABCH1 sit in the contig-wide distribution? ─────────────────
ABCH1 = (24_825_285, 24_848_323)
gene = (df.position >= ABCH1[0]) & (df.position <= ABCH1[1])
print(f"\nABCH1 windows: {gene.sum()}  (percentiles are within "
      f"{CONTIG}:{REGION_START:,}-{REGION_STOP:,}, not genome-wide)", flush=True)
print(f"{'cohort':38s} {'max PBS':>9s} {'pctile':>7s} {'median':>8s}", flush=True)
for cohort_id, _, _ in COHORTS:
    v = df[cohort_id].values
    peak = np.nanmax(v[gene.values])
    pct = 100 * np.nanmean(v[~np.isnan(v)] < peak)
    print(f"{cohort_id:38s} {peak:9.3f} {pct:7.2f} {np.nanmedian(v):8.3f}", flush=True)

print("\ndone", flush=True)
