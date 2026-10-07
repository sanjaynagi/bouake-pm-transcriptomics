---
name: malariagen-data-phenotypes-igv
description: API reference for malariagen_data phenotype data (Ag3 only) — phenotype_data(), phenotypes_with_snps(), phenotypes_with_haplotypes(), phenotype_binary(), phenotype_sample_sets(); and IGV browser integration — igv(), view_alignments(). Use when working with insecticide resistance bioassay data or visualising read alignments.
---

# malariagen_data — Phenotype data and IGV

---

## Phenotype data (Ag3 only)

Phenotype data includes results from insecticide bioassays (e.g. WHO tube tests, CDC bottle assays) linked to sequenced samples.

### `phenotype_sample_sets`

Property listing sample sets that have associated phenotype data.

```python
ag3.phenotype_sample_sets
# e.g. ['AG1000G-BF-A', 'AG1000G-GH', ...]
```

---

### `phenotype_data()`

Load all phenotype records for one or more sample sets.

```python
ag3.phenotype_data(
    sample_sets=None,            # str or list[str] — sample set IDs with phenotype data
    sample_query=None,
    sample_query_options=None,
)
```

**Returns:** `pd.DataFrame` — one row per sample × insecticide × concentration × bioassay combination. Key columns include:

| Column | Description |
|---|---|
| `sample_id` | Sample identifier |
| `insecticide` | Insecticide tested (e.g. `"permethrin"`, `"deltamethrin"`) |
| `concentration` | Insecticide concentration used |
| `phenotype` | Outcome: `"resistant"`, `"susceptible"`, or `"intermediate"` |
| `alive_dead` | Raw survival outcome |

---

### `phenotype_binary()`

Return a binary resistance/susceptibility call per sample.

```python
ag3.phenotype_binary(
    insecticide,                 # str — e.g. "permethrin"
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
)
```

**Returns:** `pd.DataFrame` with `sample_id` and `phenotype_binary` columns (1 = resistant, 0 = susceptible).

---

### `phenotypes_with_snps()`

Join phenotype data with SNP genotypes for a transcript.

```python
ag3.phenotypes_with_snps(
    transcript,                  # str — transcript ID, e.g. "AGAP004707-RA"
    insecticide,                 # str — insecticide name
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    site_mask=None,
    effects=True,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame` — one row per sample, with phenotype columns and one column per variant.

---

### `phenotypes_with_haplotypes()`

Join phenotype data with haplotype cluster assignments.

```python
ag3.phenotypes_with_haplotypes(
    region,
    insecticide,
    analysis=base_params.DEFAULT,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `pd.DataFrame`

---

## IGV — Integrative Genomics Viewer

These functions launch an embedded IGV browser in a Jupyter notebook to visualise read alignments.

### `igv()`

Launch a basic IGV browser at a given locus.

```python
ag3.igv(
    region,                      # str — genomic locus, e.g. "2L:28_000_000-28_100_000"
    tracks=None,                 # list[dict] — custom IGV track configurations
)
```

**Returns:** `igv-notebook` browser widget (displayed inline in Jupyter).

---

### `view_alignments()`

Launch an IGV browser showing read alignments for one or more samples.

```python
ag3.view_alignments(
    region,
    sample_id,                   # str or list[str] — sample(s) to display
    sample_sets=None,
    visibility_window=20_000,    # int — window size (bp) below which reads become visible
)
```

**Returns:** `igv-notebook` browser widget.

**Notes:**
- Requires access credentials to the cloud storage bucket — run inside Google Colab or with valid GCP credentials.
- `visibility_window` controls at what zoom level read-level data is loaded. Reduce for large regions.
