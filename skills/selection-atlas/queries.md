---
name: selection-atlas-queries
description: Recipes for interrogating the Selection Atlas signals table — find genes overlapping sweeps, rank cohorts by signal at a locus, find parallel selection across cohorts/species, summarise sweep strength and breadth, identify novel (non-canonical) signals. Use when asked about which cohorts have a sweep at a gene, which genes are under sweeps, how strong the sweeps are, or for any other filtering / ranking task on the signals CSV.
---

# Selection Atlas — Common queries

All examples assume:

```python
import pandas as pd
import numpy as np

signals = pd.read_csv("h12-signal-detection-all.csv")
```

For genome-feature lookups, `malariagen_data` is used (see also `skills-malariagen_data/`):

```python
import malariagen_data
ag3 = malariagen_data.Ag3()
# af1 = malariagen_data.Af1()  # for An. funestus
```

---

## 1. Which cohorts have a sweep at a given gene?

Given an AGAP gene identifier, find every cohort with a signal whose `span1` or `span2` covers the gene.

```python
gene_id = "AGAP004707"  # Vgsc

gene = ag3.geneset().query("ID == @gene_id").iloc[0]
contig = gene["contig"]
gstart, gstop = int(gene["start"]), int(gene["end"])

hits = signals.query(
    "contig == @contig "
    "and span1_pstart <= @gstop and span1_pstop >= @gstart"
).sort_values("delta_i", ascending=False)

hits[["cohort_id", "pcenter", "stat_max", "delta_i", "decay"]]
```

Change `span1_*` → `focus_*` for the tightest association, or `span2_*` for the broadest.

---

## 2. Which genes fall inside a given signal?

Given a single signal (one row), return all protein-coding genes overlapping its `span1` interval.

```python
sig = signals.iloc[0]
genes = ag3.geneset(attributes=["ID", "Name", "description"]).query(
    "contig == @sig.contig and type == 'gene' "
    "and start <= @sig.span1_pstop and end >= @sig.span1_pstart"
)
genes[["ID", "Name", "description", "start", "end"]]
```

For a whole cohort's signals, iterate / merge with `pd.merge_asof` on interval bounds, or use `pyranges`.

---

## 3. Top N strongest signals overall

Rank by **confidence** (`delta_i`) or **intensity** (`stat_max`).

```python
# most confident peaks (most peak-like vs flat)
signals.nlargest(20, "delta_i")[
    ["cohort_id", "contig", "pcenter", "stat_max", "delta_i"]
]

# most intense sweeps (closest to a complete sweep)
signals.nlargest(20, "stat_max")[
    ["cohort_id", "contig", "pcenter", "stat_max", "delta_i"]
]
```

---

## 4. Signals unique to / shared across taxa

The taxon is encoded in `cohort_id`. Extract it:

```python
signals["taxon"] = signals["cohort_id"].str.extract(r"_(gamb|colu|arab|biss)(?:_|$)")
signals["country"] = signals["cohort_id"].str.slice(0, 2)
signals["year"] = signals["cohort_id"].str.extract(r"_(\d{4})(?:_|$)").astype(int)
```

Cluster signals at the same locus across cohorts by binning on `pcenter` (e.g. 1 Mb bins) or by overlap of `span1`:

```python
# crude 1 Mb binning
signals["locus_bin"] = (signals["pcenter"] // 1_000_000).astype(int)
parallel = (
    signals.groupby(["contig", "locus_bin"])["taxon"]
    .nunique()
    .sort_values(ascending=False)
)
parallel.head(20)  # loci under selection in the most taxa
```

For a proper overlap-based clustering, sort by `contig` then `pcenter` and merge rows whose `span1` intervals touch.

---

## 5. Parallel selection at a named locus

"Who else has a sweep at *Cyp6p* / *Gste* / *Ace1* / *Keap1* ...?"

```python
def cohorts_at_locus(contig, start, stop, interval="span1"):
    a, b = f"{interval}_pstart", f"{interval}_pstop"
    return signals.query(
        "contig == @contig and " f"{a} <= @stop and {b} >= @start"
    ).sort_values("delta_i", ascending=False)

# Example: Cyp6p cluster on 2RL
cohorts_at_locus("2RL", 28_480_000, 28_600_000)
```

---

## 6. Strongest cohort per locus

Collapse a locus cluster to one "champion" cohort:

```python
champ = (
    cohorts_at_locus("2RL", 28_480_000, 28_600_000)
    .sort_values("delta_i", ascending=False)
    .drop_duplicates(subset=["contig"], keep="first")
)
```

---

## 7. Sweep breadth and skew

Peak breadth in bp (from `span2`), and which side it decays more slowly on:

```python
signals["breadth_bp"]  = signals["span2_pstop"] - signals["span2_pstart"]
signals["asymmetry"]   = signals["decay_right"] - signals["decay_left"]   # +ve = slower right decay
signals.groupby("contig")["breadth_bp"].describe()
```

A strong positive `asymmetry` can hint at the selected variant sitting toward the **left** flank of the signal (the longer decay is away from the variant due to recombination differences).

---

## 8. Novel vs canonical signals

Canonical insecticide-resistance loci (from the manuscript) that you may want to exclude or flag:

| Name | Contig | Approx. region (bp, *Ag3*) |
|---|---|---|
| Vgsc  | 2RL | ~2,358,000–2,430,000  (AGAP004707) |
| Rdl   | 2RL | ~25,360,000–25,440,000 (AGAP006028) |
| Cyp6p cluster | 2RL | ~28,480,000–28,600,000 |
| Gste cluster  | 3RL | ~28,590,000–28,620,000 |
| Ace1   | 3RL | ~3,490,000–3,510,000 (AGAP001356) |
| Cyp9k1 | X   | ~15,240,000–15,250,000 (AGAP000818) |
| Coeaexf| 2RL | Carboxylesterase (novel) |
| Keap1  | 3RL | Metabolic regulator (novel in *gambiae* s.l.) |

Tag signals with the best-matching canonical locus (or `"novel"`):

```python
known = [
    ("Vgsc",   "2RL",  2_350_000,  2_440_000),
    ("Rdl",    "2RL", 25_340_000, 25_450_000),
    ("Cyp6p",  "2RL", 28_400_000, 28_700_000),
    ("Gste",   "3RL", 28_500_000, 28_700_000),
    ("Ace1",   "3RL",  3_450_000,  3_550_000),
    ("Cyp9k1", "X",   15_200_000, 15_300_000),
]

def label(row):
    for name, c, a, b in known:
        if row.contig == c and row.span1_pstart <= b and row.span1_pstop >= a:
            return name
    return "novel"

signals["locus"] = signals.apply(label, axis=1)
signals.groupby("locus").size().sort_values(ascending=False)
```

Verify exact coordinates by looking up the gene in `ag3.geneset()` rather than trusting these approximate ranges.

---

## 9. Aggregating across cohorts

Count distinct cohorts per locus, per country, per year:

```python
agg = (
    signals.groupby("locus")
    .agg(
        n_cohorts=("cohort_id", "nunique"),
        n_countries=("country", "nunique"),
        median_stat_max=("stat_max", "median"),
        median_delta_i=("delta_i", "median"),
        min_year=("year", "min"),
        max_year=("year", "max"),
    )
    .sort_values("n_cohorts", ascending=False)
)
```

---

## 10. Temporal trends at a locus

Across years within one country/taxon — e.g. is a sweep growing in intensity?

```python
cyp6p = cohorts_at_locus("2RL", 28_480_000, 28_600_000).copy()
cyp6p["year"] = cyp6p["cohort_id"].str.extract(r"_(\d{4})(?:_|$)").astype(int)
cyp6p.groupby("year")["stat_max"].mean()
```
