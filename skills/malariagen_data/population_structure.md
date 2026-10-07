---
name: malariagen-data-population-structure
description: API reference for malariagen_data population structure functions — pca(), plot_pca_*(), njt(), plot_njt(), biallelic_diplotype_pairwise_distances(), fst_gwss(), average_fst(), pairwise_average_fst(), plot_pairwise_average_fst(), plot_fst_gwss(), plot_diplotype_clustering(), plot_haplotype_clustering(), plot_haplotype_network(). Use when analysing genetic differentiation, clustering, or visualising population relationships.
---

# malariagen_data — Population structure

---

## Principal components analysis (PCA)

### `pca()`

```python
ag3.pca(
    region,                              # str or list[str] — genomic region(s)
    n_snps,                              # int — target number of SNPs (after thinning)
    n_components=20,                     # int — number of PCs to compute
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    min_minor_ac=2,                      # int — min minor allele count to retain a site
    max_missing_an=0,                    # int — max missing allele number per site
    imputation_method="mean",            # str — how to impute missing genotypes
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    exclude_samples=None,                # str or list[str] — sample IDs to exclude from all samples
    fit_exclude_samples=None,            # str or list[str] — exclude only from PCA fitting (but project them)
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(df_pca, evr)`
- `df_pca` — `pd.DataFrame` of shape `(samples, n_components)`. Index is `sample_id`. Column names are `PC1`, `PC2`, ...
- `evr` — `np.ndarray` of explained variance ratios, one per component.

**Notes:**
- `n_snps` controls random thinning to reduce linkage between variants. Actual SNP count after thinning may differ.
- `fit_exclude_samples` is useful for projecting admixed or related individuals onto a PCA fit without those individuals influencing the axes.

**Example:**
```python
df_pca, evr = ag3.pca(
    region="3L",
    n_snps=100_000,
    sample_sets="3.0",
    sample_query="taxon in ('gambiae', 'coluzzii')",
)
```

---

### `plot_pca_variance()`

Scree plot of explained variance.

```python
ag3.plot_pca_variance(
    evr,                    # np.ndarray — from pca()
    n_components=8,
    show=True,
)
```

**Returns:** Plotly figure.

---

### `plot_pca_coords()`

2-D scatter plot of PCA coordinates.

```python
ag3.plot_pca_coords(
    df_pca,                 # pd.DataFrame — from pca()
    x="PC1",
    y="PC2",
    color=None,             # str — metadata column to colour points
    symbol=None,            # str — metadata column to vary point symbol
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    title=None,
    show=True,
)
```

**Returns:** Plotly figure.

---

### `plot_pca_coords_3d()`

3-D scatter plot of PCA coordinates.

```python
ag3.plot_pca_coords_3d(
    df_pca,
    x="PC1",
    y="PC2",
    z="PC3",
    color=None,
    symbol=None,
    # ... same as plot_pca_coords()
)
```

**Returns:** Plotly figure.

---

## Genetic distance and neighbour-joining trees

### `biallelic_diplotype_pairwise_distances()`

Compute pairwise genetic distances between samples.

```python
ag3.biallelic_diplotype_pairwise_distances(
    region,
    n_snps,
    metric="cityblock",          # str — distance metric (e.g. "cityblock", "euclidean")
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    min_minor_ac=2,
    max_missing_an=0,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(dist, samples)` — `np.ndarray` distance matrix of shape `(n, n)` and list of sample IDs.

---

### `njt()`

Compute a neighbour-joining tree from a distance matrix.

```python
ag3.njt(
    region,
    n_snps,
    # ... same params as biallelic_diplotype_pairwise_distances()
)
```

**Returns:** `(tree, samples)` — scikit-bio `TreeNode` and sample IDs.

---

### `plot_njt()`

Plot the neighbour-joining tree.

```python
ag3.plot_njt(
    region,
    n_snps,
    # ... same params
    color=None,             # str — metadata column to colour leaves
    symbol=None,
    sample_sets=None,
    sample_query=None,
    title=None,
    show=True,
)
```

**Returns:** Plotly figure.

---

## Fst analysis

### `fst_gwss()`

Genome-wide Fst scan between two cohorts.

```python
ag3.fst_gwss(
    contig,
    window_size,                 # int — number of SNPs per window
    cohort1_query,               # str — pandas query for cohort 1
    cohort2_query,               # str — pandas query for cohort 2
    sample_query_options=None,
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    cohort_size=None,
    min_cohort_size=15,
    max_cohort_size=50,
    random_seed=42,
    inline_array=True,
    chunks="native",
    clip_min=0.0,                # float — clip negative Fst values to this minimum
)
```

**Returns:** `(x, fst)` — position and Fst arrays.

### `plot_fst_gwss()`

```python
ag3.plot_fst_gwss(
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

### `average_fst()`

Compute genome-wide average Fst between two cohorts.

```python
ag3.average_fst(
    region,
    cohort1_query,
    cohort2_query,
    sample_query_options=None,
    sample_sets=None,
    site_mask=base_params.DEFAULT,
    cohort_size=None,
    min_cohort_size=15,
    max_cohort_size=50,
    random_seed=42,
    n_jack=200,                  # int — number of jackknife blocks for SE estimation
    confidence_level=0.95,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` with columns `fst`, `se`, `ci_lower`, `ci_upper`.

---

### `pairwise_average_fst()`

Compute average Fst for all pairs of cohorts.

```python
ag3.pairwise_average_fst(
    region,
    cohorts,                     # str or dict — cohort grouping
    min_cohort_size=10,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    site_mask=base_params.DEFAULT,
    n_jack=200,
    confidence_level=0.95,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — pairwise Fst matrix.

### `plot_pairwise_average_fst()`

```python
ag3.plot_pairwise_average_fst(
    df_fst,                      # pd.DataFrame — from pairwise_average_fst()
    title=None,
    show=True,
)
```

**Returns:** Plotly figure.

---

## Diplotype clustering

### `plot_diplotype_clustering()`

Hierarchical clustering of diplotypes (samples) based on genotype distance.

```python
ag3.plot_diplotype_clustering(
    region,
    n_snps,
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    color=None,                  # str — metadata column for leaf colouring
    linkage_method="complete",
    count_sort=True,
    distance_sort=False,
    title=None,
    show=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** Plotly figure.

### `plot_diplotype_clustering_advanced()`

Advanced version with cohort-level grouping and metadata overlays.

```python
ag3.plot_diplotype_clustering_advanced(
    region,
    n_snps,
    cohort1_query,
    cohort2_query=None,
    # ... same params
)
```

---

## Haplotype clustering and network analysis

### `plot_haplotype_clustering()`

Hierarchical clustering of haplotypes.

```python
ag3.plot_haplotype_clustering(
    region,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    color=None,
    linkage_method="complete",
    count_sort=True,
    distance_sort=False,
    title=None,
    show=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** Plotly figure.

---

### `haplotype_pairwise_distances()`

Compute pairwise distances between haplotypes.

```python
ag3.haplotype_pairwise_distances(
    region,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `(dist, haplotypes)` — distance matrix and haplotype labels.

---

### `plot_haplotype_network()`

Median-joining haplotype network visualisation.

```python
ag3.plot_haplotype_network(
    region,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    color=None,
    max_dist=2,                  # int — max pairwise distance for connecting haplotypes
    show=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** Bokeh figure.
