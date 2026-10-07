---
name: malariagen-data-snp-data
description: API reference for malariagen_data SNP data access — snp_calls(), snp_allele_counts(), site_annotations(), is_accessible(), biallelic_snp_calls(), biallelic_diplotypes(), biallelic_snps_to_plink(), plot_snps(). Use when loading genotype data or preparing SNP datasets for analysis.
---

# malariagen_data — SNP data

---

## `site_mask_ids`

Property listing available site mask identifiers.

```python
ag3.site_mask_ids
# e.g. ('gamb_colu', 'gamb_colu_arab', 'arab')
```

---

## `snp_calls()`

Access SNP genotype calls for one or more genomic regions.

```python
ag3.snp_calls(
    region,                     # str or list[str] — contig(s) or region(s)
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=None,             # str — e.g. "gamb_colu_arab"; None = no filter
    site_class=None,            # str — e.g. "CDS_DEG_4" for 4-fold degenerate sites
    inline_array=True,
    chunks="native",
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
)
```

**Returns:** `xarray.Dataset` with dimensions `(variants, samples)`. Key variables:

| Variable | Shape | Description |
|---|---|---|
| `variant_position` | `(variants,)` | 1-based genomic position |
| `variant_allele` | `(variants, 4)` | REF + up to 3 ALT alleles |
| `variant_contig` | `(variants,)` | Contig index |
| `call_genotype` | `(variants, samples, 2)` | Diploid genotype calls (0=REF, 1-3=ALT, -1=missing) |
| `call_genotype_phased` | `(variants, samples)` | Whether calls are phased |
| `call_GQ` | `(variants, samples)` | Genotype quality |
| `call_AD` | `(variants, samples, 4)` | Allele depth |
| `sample_id` | `(samples,)` | Sample identifiers |

**Example:**
```python
ds = ag3.snp_calls(
    region="2L:28_000_000-29_000_000",
    sample_sets="3.0",
    sample_query="country == 'Ghana'",
    site_mask="gamb_colu_arab",
)
```

---

## `snp_allele_counts()`

Count allele occurrences across samples for a genomic region.

```python
ag3.snp_allele_counts(
    region,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=None,
    site_class=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `np.ndarray` of shape `(variants, 4)` — allele counts per site (REF, ALT1, ALT2, ALT3).

---

## `site_annotations()`

Load functional site annotations (codon position, amino acid effect, etc.).

```python
ag3.site_annotations(
    region,
    site_mask=None,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` with variables including `variant_contig`, `variant_position`, `codon_position`, `codon_degeneracy`, `codon_nonsyn`, `codon_stop`, `seq_cls`, `seq_flen`, `seq_relpos_start`, `seq_relpos_stop`.

---

## `is_accessible()`

Boolean array indicating whether each site passes the site mask.

```python
ag3.is_accessible(
    region,                # str — contig or region
    site_mask,             # str — required; e.g. "gamb_colu_arab"
    inline_array=True,
    chunks="native",
)
```

**Returns:** `np.ndarray` of `bool`, length = number of sites in region.

---

## `biallelic_snp_calls()`

Filter SNP calls to biallelic sites only (exactly one ALT allele observed), optionally LD-pruned.

```python
ag3.biallelic_snp_calls(
    region,
    n_snps,                    # int — target number of SNPs after thinning
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    min_minor_ac=None,         # int — minimum minor allele count to retain a site
    max_missing_an=None,       # int — max number of missing alleles allowed per site
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` — same structure as `snp_calls()` but restricted to biallelic sites.

---

## `biallelic_diplotypes()`

Return a 2-D array of diplotype integers (allele dosage 0/1/2) for biallelic sites.

```python
ag3.biallelic_diplotypes(
    region,
    n_snps,
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    min_minor_ac=None,
    max_missing_an=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `np.ndarray` of shape `(samples, variants)` with values 0, 1, or 2.

---

## `biallelic_snps_to_plink()`

Export biallelic SNPs to PLINK BED/BIM/FAM format.

```python
ag3.biallelic_snps_to_plink(
    region,
    n_snps,
    output_prefix,             # str — path prefix for output files
    thin_offset=0,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    site_class=None,
    min_minor_ac=None,
    max_missing_an=None,
    cohort_size=None,
    min_cohort_size=None,
    max_cohort_size=None,
    random_seed=42,
    inline_array=True,
    chunks="native",
)
```

Writes `.bed`, `.bim`, `.fam` files. No return value.

---

## `plot_snps()`

Plot SNP density along a genomic region.

```python
ag3.plot_snps(
    region,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    site_mask=base_params.DEFAULT,
    sizing_mode="stretch_width",
    width=800,
    show=True,
    output_backend="svg",
    x_range=None,
)
```

**Returns:** Bokeh figure.
