---
name: anoexpress-data
description: API reference for AnoExpress data-loading functions — ae.data(), ae.metadata(), ae.sample_metadata(), ae.xpress_metadata(), ae.irtex_metadata(), ae.load_annotations(), ae.load_results_arrays(). Use when loading expression data or metadata from the AnoExpress meta-analysis.
---

# AnoExpress — Data loading

All functions are importable from the top-level `anoexpress` namespace.

---

## `ae.data()`

The primary entry point for loading expression data.

```python
ae.data(
    data_type,           # "fcs" | "pvals" | "log2counts"
    analysis,            # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun"
    microarray=False,    # include IR-Tex microarray data
    gene_id=None,        # str, list, or file path
    sample_query=None,   # pandas query string to subset comparisons
    sort_by=None,        # "median" | "mean" | "agap" | "position" | None
    annotations=False,   # add GeneName/GeneDescription to index
    pvalue_filter=None,  # float — null FCs where padj > threshold (fcs only)
    low_count_filter=None,  # int — remove genes with median counts below threshold
    fraction_na_allowed=None,  # float 0–1, max fraction of NaN per gene
    gff_method='vectorbase',   # "vectorbase" | "malariagen_data"
)
```

**Returns:** `pd.DataFrame` — rows are genes (GeneID index), columns are comparison IDs (fcs/pvals) or sample IDs (log2counts).

**Notes:**
- `pvalue_filter` is only applied when `data_type="fcs"`. It sets fold-change entries to `NaN` where the corresponding adjusted p-value exceeds the threshold.
- `sort_by="position"` requires `analysis != "fun"` and downloads a GFF file (cached at `~/.cache/anoexpress/`).
- When `annotations=True` the index becomes a MultiIndex of `(GeneID, GeneName, GeneDescription)`.

**Example:**
```python
import anoexpress as ae

# fold changes for two genes, filtered by significance
fc = ae.data(
    data_type="fcs",
    analysis="gamb_colu_arab_fun",
    gene_id=["AGAP004707", "AGAP002865"],
    pvalue_filter=0.05,
    annotations=True,
)
```

---

## `ae.metadata()`

Load comparison-level metadata (one row per resistant vs susceptible comparison).

```python
ae.metadata(
    analysis,          # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun"
    microarray=False,  # whether to include IR-Tex microarray comparisons
)
```

**Returns:** `pd.DataFrame` with columns including `comparison`, `species`, `country`, `resistant`, `susceptible`, `technology`.

---

## `ae.sample_metadata()`

Load sample-level metadata (one row per RNA-Seq sample).

```python
ae.sample_metadata(analysis)
```

**Returns:** `pd.DataFrame` with columns including `sampleID`, `species`, `country`, `condition`.

---

## `ae.xpress_metadata()`

Load the raw AnoExpress comparison metadata (all species, unfiltered).

```python
ae.xpress_metadata()
```

**Returns:** `pd.DataFrame`

---

## `ae.irtex_metadata()`

Load IR-Tex microarray comparison metadata.

```python
ae.irtex_metadata()
```

**Returns:** `pd.DataFrame`

---

## `ae.load_annotations()`

Load PFAM and GO annotations for *Anopheles gambiae* proteins.

```python
ae.load_annotations()
```

**Returns:** `pd.DataFrame` with columns `transcript`, `pstart`, `pend`, `pfamid`, `domain`, `domseq`, `GO_terms`, `gene_id`.

---

## `ae.load_results_arrays()`

Low-level loader for a single data type / analysis combination directly from GitHub.

```python
ae.load_results_arrays(
    data_type,   # "fcs" | "pvals" | "log2counts"
    analysis,    # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun" | "irtex"
)
```

**Returns:** `pd.DataFrame` indexed by `GeneID`.

---

## Internal helpers (not part of public API)

| Function | Purpose |
|---|---|
| `filter_low_counts()` | Removes genes below a median count threshold |
| `filter_nas()` | Removes genes exceeding a missing-value fraction |
| `null_fold_changes()` | Sets FC values to NaN where padj > threshold |
| `add_annotations_to_array()` | Merges gene name/description into a data frame |
| `_sort_genes()` | Sorts rows by median, mean, AGAP ID, or genomic position |
