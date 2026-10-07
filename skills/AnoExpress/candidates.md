---
name: anoexpress-candidates
description: API reference for AnoExpress candidate gene functions — ae.load_candidates(), ae.consistent_genes(), ae.contig_expression(). Use when ranking candidate genes by expression or finding consistently differentially expressed genes.
---

# AnoExpress — Candidate genes

Functions for identifying and ranking candidate insecticide resistance genes. Importable from the top-level `anoexpress` namespace.

---

## `ae.load_candidates()`

Rank all genes by a summary statistic (e.g. median fold change across experiments) and return a sorted dataframe. Optionally filter by annotation or fold change threshold.

```python
ae.load_candidates(
    analysis,               # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun"
    name='median',          # label for the ranking column in the output
    func=np.nanmedian,      # ranking function — any numpy-compatible row aggregator
    query_annotation=None,  # GO term or PFAM domain to pre-filter genes
    query_fc=None,          # float — keep only genes with fold change > this value
    microarray=False,
    low_count_filter=None,  # int — remove genes with median counts below threshold
    fraction_na_allowed=None,  # float 0–1
)
```

**Returns:** `pd.DataFrame` with columns:
- `GeneID`, `GeneName`, `GeneDescription`
- `{name} log2 Fold Change` — raw log2 value
- `{name} Fold Change` — back-transformed (rounded to 2 d.p.)

Rows are sorted descending by the summary statistic.

**Example:**
```python
import numpy as np
import anoexpress as ae

# top candidates by median fold change, keeping only those > 2-fold up
candidates = ae.load_candidates(
    analysis="gamb_colu_arab_fun",
    func=np.nanmedian,
    name="median",
    query_fc=2,
)
print(candidates.head(20))
```

**Example — filter by PFAM domain:**
```python
gst_candidates = ae.load_candidates(
    analysis="gamb_colu_arab_fun",
    query_annotation="GST_N",
    func=np.nanmedian,
)
```

---

## `ae.consistent_genes()`

Find genes that are significantly differentially expressed in at least `n_experiments` comparisons.

```python
ae.consistent_genes(
    analysis,        # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun"
    direction,       # "up" | "down"
    n_experiments,   # int — minimum number of experiments required
    low_count_filter=None,
)
```

**Returns:** `pd.DataFrame` (fold change data for qualifying genes), or `None` if no genes pass the threshold.

**Logic:**
- `direction="up"`: genes where FC > 0 in ≥ `n_experiments` comparisons **and** padj < 0.05 in > `n_experiments` comparisons.
- `direction="down"`: same logic with FC < 0.

**Example:**
```python
# genes up-regulated in at least 5 experiments
up_genes = ae.consistent_genes(
    analysis="gamb_colu_arab_fun",
    direction="up",
    n_experiments=5,
)
```

---

## `ae.contig_expression()`

Calculate per-gene expression values along a chromosome arm and compute a sliding-window median. Used internally by `plot_contig_expression_track()` but also useful standalone.

```python
ae.contig_expression(
    contig,                  # "2L" | "2R" | "3L" | "3R" | "X" | "2RL" | "3RL"
    analysis,
    data_type='fcs',         # "fcs" | "log2counts"
    microarray=False,
    pvalue_filter=None,
    size=10,                 # window size in genes
    step=5,                  # step size in genes
    fraction_na_allowed=None,
)
```

**Returns:** `(fold_change_df, windowed_fold_change_df)`

- `fold_change_df` — long-format DataFrame with columns `midpoint`, `GeneID`, `GeneName`, `GeneDescription`, `comparison`, `fold_change`
- `windowed_fold_change_df` — DataFrame with columns `midpoint`, `median_fc` (sliding window statistics)

**Note:** Requires `malariagen_data` and `scikit-allel` (`allel`) to be installed.
