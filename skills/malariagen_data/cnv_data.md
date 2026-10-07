---
name: malariagen-data-cnv
description: API reference for malariagen_data CNV data — cnv_hmm(), cnv_coverage_calls(), cnv_discordant_read_calls(), gene_cnv(), plot_cnv_hmm_coverage(), plot_cnv_hmm_heatmap(). Use when working with copy number variation, coverage-based CNV calls, or gene amplification/deletion analysis.
---

# malariagen_data — CNV data

---

## `coverage_calls_analysis_ids`

Property listing available CNV coverage calls analysis identifiers.

```python
ag3.coverage_calls_analysis_ids
# e.g. ('gamb_colu',)
```

---

## `cnv_hmm()`

Access Hidden Markov Model (HMM)-based CNV calls based on normalised read coverage.

```python
ag3.cnv_hmm(
    region,                               # str or list[str] — contig(s) or region(s)
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    max_coverage_variance=0.2,            # float — exclude samples with high coverage variance
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` with dimensions `(variants, samples)`. Key variables:

| Variable | Shape | Description |
|---|---|---|
| `variant_position` | `(variants,)` | Start of HMM window (1-based) |
| `variant_end` | `(variants,)` | End of HMM window |
| `call_CN` | `(variants, samples)` | Called copy number (integer) |
| `call_RawCov` | `(variants, samples)` | Raw normalised coverage |
| `call_NormCov` | `(variants, samples)` | Normalised coverage |
| `sample_id` | `(samples,)` | Sample identifiers |
| `sample_coverage_variance` | `(samples,)` | Per-sample coverage variance |
| `sample_is_high_variance` | `(samples,)` | Whether sample was flagged as high-variance |

**Notes:**
- Windows are ~300 bp each.
- `max_coverage_variance` filters out samples with unstable coverage profiles that would confound CNV calling.

---

## `cnv_coverage_calls()`

Access pre-called CNV regions from a coverage-based analysis (larger, validated CNV events).

```python
ag3.cnv_coverage_calls(
    region,
    sample_set,                             # str — a single sample set ID
    analysis=base_params.DEFAULT,           # str — coverage calls analysis version
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` — CNV call intervals with copy number and confidence.

---

## `cnv_discordant_read_calls()` (Ag3 only)

Access CNV calls based on discordant read pairs (detects breakpoints).

```python
ag3.cnv_discordant_read_calls(
    contig,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset`

---

## `gene_cnv()`

Compute modal copy number per gene using HMM data.

```python
ag3.gene_cnv(
    region,                              # str or list[str] — gene region(s) of interest
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    max_coverage_variance=0.2,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` with dimensions `(genes, samples)`. Key variables:

| Variable | Description |
|---|---|
| `gene_id` | Gene identifier |
| `gene_name` | Gene name |
| `gene_contig` / `gene_start` / `gene_end` | Genomic coordinates |
| `CN_mode` | Modal copy number per gene per sample |
| `CN_mode_count` | Number of HMM windows supporting the modal CN |
| `sample_coverage_variance` | Per-sample coverage variance |

**Example:**
```python
ds = ag3.gene_cnv(
    region="2L:28_500_000-29_000_000",
    sample_sets="3.0",
    sample_query="country == 'Ghana'",
)
import pandas as pd
df = ds[["gene_id", "gene_name", "CN_mode"]].to_dataframe()
```

---

## `plot_cnv_hmm_coverage()`

Multi-panel plot of HMM copy number and raw coverage for a genomic region.

```python
ag3.plot_cnv_hmm_coverage(
    region,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    max_coverage_variance=0.2,
    sizing_mode="stretch_width",
    width=800,
    show=True,
    output_backend="svg",
)
```

**Returns:** Bokeh layout (column of figures).

---

## `plot_cnv_hmm_heatmap()`

Heatmap of CNV copy number (HMM) across samples and genomic windows.

```python
ag3.plot_cnv_hmm_heatmap(
    region,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    max_coverage_variance=0.2,
    sizing_mode="stretch_width",
    width=800,
    row_height=3,
    show=True,
    output_backend="svg",
)
```

**Returns:** Bokeh figure.
