---
name: malariagen-data-frequencies
description: API reference for malariagen_data frequency analysis — snp_allele_frequencies(), aa_allele_frequencies(), snp_allele_frequencies_advanced(), aa_allele_frequencies_advanced(), gene_cnv_frequencies(), gene_cnv_frequencies_advanced(), haplotypes_frequencies(), haplotypes_frequencies_advanced(), and plotting functions. Use when computing or visualising allele/CNV/haplotype frequencies across cohorts.
---

# malariagen_data — Frequency analysis

All frequency functions accept a `cohorts` parameter which defines how samples are grouped. It can be:
- A `str` naming a metadata column (e.g. `"admin1_year"`) — one cohort per unique value
- A `dict` mapping cohort labels to pandas query strings (e.g. `{"Ghana": "country == 'Ghana'"}`)

---

## `snp_allele_frequencies()`

Compute per-cohort SNP allele frequencies for a transcript or region.

```python
ag3.snp_allele_frequencies(
    transcript=None,             # str — transcript ID (e.g. "AGAP004707-RA"); use this OR region
    region=None,                 # str — genomic region; use this OR transcript
    cohorts=None,                # str or dict — cohort grouping
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    site_mask=None,
    sample_sets=None,
    drop_invariant=True,         # bool — remove sites with no variation across cohorts
    effects=True,                # bool — include codon effect annotations
    include_counts=False,        # bool — include raw allele counts alongside frequencies
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per variant × allele combination. Columns include `contig`, `position`, `ref_allele`, `alt_allele`, `effect`, `impact`, and one `frq_<cohort>` column per cohort.

---

## `aa_allele_frequencies()`

Compute amino-acid-level allele frequencies (SNPs grouped by protein change).

```python
ag3.aa_allele_frequencies(
    transcript,                  # str — transcript ID (required)
    cohorts=None,
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    site_mask=None,
    sample_sets=None,
    drop_invariant=True,
    include_counts=False,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per amino acid change. Includes `aa_change`, `ref_aa`, `alt_aa`, and `frq_<cohort>` columns.

---

## `snp_allele_frequencies_advanced()`

Advanced SNP frequency analysis with automatic cohort grouping by taxon × area × time.

```python
ag3.snp_allele_frequencies_advanced(
    transcript=None,
    region=None,
    area_by="admin1_iso",        # str — spatial grouping column
    period_by="year",            # str — temporal grouping column
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    site_mask=None,
    drop_invariant=True,
    effects=True,
    include_counts=False,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` — dimensions `(cohorts, variants)`. Variables prefixed with `cohort_` describe each cohort; variables prefixed with `variant_` describe each site; variables prefixed with `event_` contain allele counts and frequencies.

---

## `aa_allele_frequencies_advanced()`

Advanced amino-acid-level frequency analysis.

```python
ag3.aa_allele_frequencies_advanced(
    transcript,
    area_by="admin1_iso",
    period_by="year",
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    site_mask=None,
    drop_invariant=True,
    include_counts=False,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset`

---

## `gene_cnv_frequencies()`

Compute per-cohort frequencies of gene amplification and deletion.

```python
ag3.gene_cnv_frequencies(
    region,
    cohorts=None,
    sample_query=None,
    sample_query_options=None,
    sample_sets=None,
    min_cohort_size=10,
    max_coverage_variance=0.2,
    drop_invariant=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per gene × CNV type (amp/del). Includes `gene_id`, `gene_name`, `cnv_type`, and `frq_<cohort>` columns.

---

## `gene_cnv_frequencies_advanced()`

Advanced gene CNV frequency analysis with automatic cohort grouping.

```python
ag3.gene_cnv_frequencies_advanced(
    region,
    area_by="admin1_iso",
    period_by="year",
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    max_coverage_variance=0.2,
    drop_invariant=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset`

---

## `haplotypes_frequencies()`

Compute per-cohort haplotype frequencies for a genomic region.

```python
ag3.haplotypes_frequencies(
    region,
    cohorts,                     # str or dict — required
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    sample_sets=None,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per distinct haplotype, sorted by `max_af`. Columns are `frq_<cohort>` and `max_af`. Index is a haplotype label (`"H0"`, `"H1"`, ...).

---

## `haplotypes_frequencies_advanced()`

Advanced haplotype frequency analysis with automatic cohort grouping.

```python
ag3.haplotypes_frequencies_advanced(
    region,
    area_by="admin1_iso",
    period_by="year",
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    min_cohort_size=10,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset`

---

## Plotting frequency results

### `plot_frequencies_heatmap()`

```python
ag3.plot_frequencies_heatmap(
    df,                          # pd.DataFrame — output of a *_frequencies() function
    index=None,                  # str or list[str] — columns to use as row labels
    max_len=100,                 # int — max rows to display
    col_width=40,                # int — column width in pixels
    row_height=20,               # int — row height in pixels
    x_label="Cohort",
    y_label="Variant",
    colorscale="Reds",
    zmin=0,
    zmax=1,
    title=None,
    width=None,
    height=None,
    show=True,
    renderer=None,
)
```

**Returns:** Plotly figure.

---

### `plot_frequencies_time_series()`

```python
ag3.plot_frequencies_time_series(
    ds,                          # xarray.Dataset — output of a *_frequencies_advanced() function
    height=200,
    width=None,
    facet_col_wrap=6,
    genes=None,                  # list[str] — gene IDs to annotate (for CNV data)
    title=None,
    show=True,
    renderer=None,
)
```

**Returns:** Plotly figure.

---

### `plot_frequencies_interactive_map()`

```python
ag3.plot_frequencies_interactive_map(
    ds,                          # xarray.Dataset — output of a *_frequencies_advanced() function
    variant,                     # int or str — variant/gene index or label to display
    period=0,                    # int — time period index
    colorscale="Reds",
    basemap="street",
    width=600,
    height=500,
    show=True,
)
```

**Returns:** Plotly figure.
