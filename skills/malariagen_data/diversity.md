---
name: malariagen-data-diversity
description: API reference for malariagen_data diversity and heterozygosity functions — diversity_stats(), cohort_diversity_stats(), plot_diversity_stats(), plot_heterozygosity(), roh_hmm(), plot_roh(). Use when computing nucleotide diversity, heterozygosity, or runs of homozygosity (ROH).
---

# malariagen_data — Diversity and heterozygosity

---

## Diversity statistics

### `diversity_stats()`

Compute nucleotide diversity (π), Watterson's θ, and Tajima's D for a region.

```python
ag3.diversity_stats(
    region,                              # str or list[str] — genomic region(s)
    cohort1_query,                       # str — pandas query string for cohort 1
    cohort2_query=None,                  # str — optional second cohort (for divergence stats)
    sample_query_options=None,
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    cohort_size=None,
    min_cohort_size=15,
    max_cohort_size=50,
    random_seed=42,
    n_jack=200,                          # int — jackknife blocks for SE estimation
    confidence_level=0.95,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` with columns including `theta_pi`, `theta_w`, `tajima_d`, and confidence interval columns (`ci_lower_*`, `ci_upper_*`). When `cohort2_query` is provided, also includes `dxy` (absolute divergence) and `fst`.

---

### `cohort_diversity_stats()`

Compute diversity statistics across multiple cohorts in one call.

```python
ag3.cohort_diversity_stats(
    cohorts,                             # str or dict — cohort grouping
    region,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    cohort_size=None,
    min_cohort_size=15,
    max_cohort_size=50,
    random_seed=42,
    n_jack=200,
    confidence_level=0.95,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per cohort, same columns as `diversity_stats()`.

---

### `plot_diversity_stats()`

Bar chart of diversity statistics across cohorts.

```python
ag3.plot_diversity_stats(
    df,                                  # pd.DataFrame — from cohort_diversity_stats()
    stat,                                # str — e.g. "theta_pi", "tajima_d"
    color=None,                          # str — metadata column for bar colours
    title=None,
    show=True,
)
```

**Returns:** Plotly figure.

---

## Heterozygosity

### `plot_heterozygosity()`

Plot per-sample heterozygosity along a genomic region.

```python
ag3.plot_heterozygosity(
    region,
    sample_id,                           # str or list[str] — sample(s) to plot
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    window_size=20_000,                  # int — window size in base pairs
    y_max=0.1,
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

## Runs of homozygosity (ROH)

### `roh_hmm()`

Detect runs of homozygosity using a Hidden Markov Model.

```python
ag3.roh_hmm(
    sample_id,                           # str or list[str] — sample(s) to analyse
    contig,                              # str — single contig
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    window_size=20_000,                  # int — window size in base pairs for pre-smoothing
    heterozygosity_height=0.15,          # float — heterozygosity level used as HMM emission parameter
    phet_roh=0.001,                      # float — prob. of heterozygosity within ROH
    phet_nonroh=(0.0025, 0.01),          # tuple — prob. of het outside ROH (multiple non-ROH states)
    transition=1e-6,                     # float — HMM state transition probability
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per detected ROH interval. Columns: `sample_id`, `contig`, `start`, `stop`, `length`, `is_roh`.

---

### `plot_roh()`

Plot ROH calls and underlying heterozygosity for one or more samples.

```python
ag3.plot_roh(
    sample_id,
    contig,
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    window_size=20_000,
    heterozygosity_height=0.15,
    phet_roh=0.001,
    phet_nonroh=(0.0025, 0.01),
    transition=1e-6,
    sizing_mode="stretch_width",
    width=800,
    height=400,
    show=True,
    output_backend="svg",
    inline_array=True,
    chunks="native",
)
```

**Returns:** Bokeh layout.
