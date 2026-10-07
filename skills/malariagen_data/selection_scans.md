---
name: malariagen-data-selection-scans
description: API reference for malariagen_data genome-wide selection scan functions — h12_gwss(), h1x_gwss(), g123_gwss(), ihs_gwss(), xpehh_gwss(), and their calibration and plotting counterparts. Use when running haplotype-based or SNP-based selection scans.
---

# malariagen_data — Genome-wide selection scans

All scan functions return raw statistics; use the corresponding `plot_*_gwss()` functions for visualisation. Results are cached to disk when `results_cache` is set at API initialisation.

---

## H12 — Extended haplotype homozygosity (single cohort)

### `h12_gwss()`

```python
ag3.h12_gwss(
    contig,                              # str — single contig, e.g. "2L"
    window_size,                         # int — number of SNPs per sliding window
    analysis=base_params.DEFAULT,        # str — phasing analysis ID
    sample_query=None,
    sample_query_options=None,
    sample_sets=None,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(x, h12, contigs)` — three `np.ndarray` arrays of genomic positions, H12 values, and contig indices.

---

### `h12_calibration()`

Calibrate the window size for H12 by estimating a null expectation.

```python
ag3.h12_calibration(
    contig,
    analysis=base_params.DEFAULT,
    sample_query=None,
    sample_query_options=None,
    sample_sets=None,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    window_sizes=(100, 200, 500, 1000, 2000),
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame`

### `plot_h12_calibration()`

```python
ag3.plot_h12_calibration(
    contig,
    analysis=base_params.DEFAULT,
    sample_query=None,
    sample_sets=None,
    # ... same params as h12_calibration()
    title=None,
    show=True,
)
```

---

### `plot_h12_gwss()`

```python
ag3.plot_h12_gwss(
    contig,
    window_size,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    title=None,
    sizing_mode="stretch_width",
    width=800,
    height=200,
    show=True,
    output_backend="svg",
    x_range=None,
)
```

**Returns:** Bokeh figure.

---

### `plot_h12_gwss_multi_panel()` / `plot_h12_gwss_multi_overlay()`

Run and plot H12 scans for multiple cohorts in separate panels or overlaid on one plot.

```python
ag3.plot_h12_gwss_multi_panel(
    contig,
    window_size,
    cohort_queries,              # list[str] — one query per cohort
    # ... shared params
)

ag3.plot_h12_gwss_multi_overlay(
    contig,
    window_size,
    cohort_queries,
    # ... shared params
)
```

---

## H1X — Extended haplotype homozygosity (between two cohorts)

### `h1x_gwss()`

```python
ag3.h1x_gwss(
    contig,
    window_size,
    cohort1_query,               # str — pandas query for cohort 1
    cohort2_query,               # str — pandas query for cohort 2
    sample_query_options=None,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(x, h1x)` — position and H1X statistic arrays.

### `plot_h1x_gwss()`

```python
ag3.plot_h1x_gwss(
    contig,
    window_size,
    cohort1_query,
    cohort2_query,
    # ... shared params
    title=None,
    show=True,
)
```

---

## G123 — Genotype-based selection scan

### `g123_gwss()`

```python
ag3.g123_gwss(
    contig,
    window_size,
    sites,                       # str — site class, e.g. "CDS_DEG_4" or a site mask ID
    site_mask=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(x, g123)` arrays.

### `g123_calibration()` / `plot_g123_calibration()`

Same pattern as H12 calibration — run across multiple window sizes to calibrate.

### `plot_g123_gwss()`

```python
ag3.plot_g123_gwss(
    contig,
    window_size,
    sites,
    # ... shared params
    title=None,
    show=True,
)
```

---

## IHS — Integrated haplotype score

### `ihs_gwss()`

```python
ag3.ihs_gwss(
    contig,
    window_size,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    site_mask=base_params.DEFAULT,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    standardize=True,            # bool — standardise IHS to Z-scores within frequency bins
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(x, ihs)` — position and |IHS| arrays.

### `plot_ihs_gwss()`

```python
ag3.plot_ihs_gwss(
    contig,
    window_size,
    # ... shared params
    title=None,
    show=True,
)
```

---

## XPEHH — Cross-population extended haplotype homozygosity

### `xpehh_gwss()`

```python
ag3.xpehh_gwss(
    contig,
    window_size,
    cohort1_query,
    cohort2_query,
    sample_query_options=None,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    cohort_size=None,
    min_cohort_size=20,
    max_cohort_size=50,
    random_seed=42,
    standardize=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(x, xpehh)` arrays.

### `plot_xpehh_gwss()`

```python
ag3.plot_xpehh_gwss(
    contig,
    window_size,
    cohort1_query,
    cohort2_query,
    # ... shared params
    title=None,
    show=True,
)
```

---

## Notes

- All GWSS functions use haplotype data (from `haplotypes()`), except `g123_gwss()` which uses diplotype SNP calls.
- Results depend on `cohort_size` subsampling — set a `random_seed` for reproducibility.
- Use `results_cache` in the API constructor to avoid re-running expensive scans.
- Typically run one contig at a time; concatenate results manually if scanning multiple chromosomes.
