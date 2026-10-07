---
name: malariagen-data-karyotype
description: API reference for malariagen_data inversion karyotyping — karyotype(). Ag3 only. Use when inferring chromosomal inversion genotypes (2La, 2Rb, etc.) from tag SNPs.
---

# malariagen_data — Inversion karyotypes (Ag3 only)

---

## `karyotype()`

Infer inversion karyotypes from tag SNPs.

```python
ag3.karyotype(
    inversion,              # str — e.g. "2La", "2Rb", "2Rc"
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
)
```

**Returns:** `pd.DataFrame` — one row per sample, with columns for homozygous/heterozygous inversion genotype.

**Available inversions:** `"2La"`, `"2Rb"`, `"2Rc"`, `"2Rd"`, `"2Rj"`, `"2Rk"` (check the Ag3 documentation for the current list).
