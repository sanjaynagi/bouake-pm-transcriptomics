---
name: malariagen-data-haplotypes
description: API reference for malariagen_data haplotype data — phasing_analysis_ids, haplotypes(), haplotype_sites(). Use when loading phased haplotype calls for a genomic region.
---

# malariagen_data — Haplotype data

---

## `phasing_analysis_ids`

Property listing available phasing analysis identifiers.

```python
ag3.phasing_analysis_ids
# e.g. ('gamb_colu', 'gamb_colu_arab', 'arab')
```

---

## `haplotypes()`

Access phased haplotype data.

```python
ag3.haplotypes(
    region,                      # str or list[str] — contig(s) or region(s)
    analysis=base_params.DEFAULT, # str — phasing analysis ID; default depends on API
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

**Returns:** `xarray.Dataset` with dimensions `(variants, samples)`. Key variables:

| Variable | Shape | Description |
|---|---|---|
| `variant_position` | `(variants,)` | 1-based position |
| `variant_allele` | `(variants, 2)` | REF and ALT alleles |
| `call_genotype` | `(variants, samples, 2)` | Phased haplotype calls (0/1 per haplotype) |
| `sample_id` | `(samples,)` | Sample identifiers |

**Notes:**
- Not all samples pass phasing QC — `haplotypes()` returns only samples with phased calls.
- The `analysis` parameter determines which phasing model/reference panel was used.

**Example:**
```python
ds_haps = ag3.haplotypes(
    region="2L:28_000_000-29_000_000",
    analysis="gamb_colu_arab",
    sample_sets="3.0",
    sample_query="country == 'Ghana'",
)
```

---

## `haplotype_sites()`

Access the list of variant sites used in a phasing analysis.

```python
ag3.haplotype_sites(
    region=None,
    analysis=base_params.DEFAULT,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset` with `variant_position`, `variant_contig`, `variant_allele`.

