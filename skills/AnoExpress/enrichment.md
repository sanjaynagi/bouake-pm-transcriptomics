---
name: anoexpress-enrichment
description: API reference for AnoExpress gene set enrichment functions — ae.go_hypergeometric(), ae.pfam_hypergeometric(), ae.kegg_hypergeometric(). Use when testing for enrichment of GO terms, PFAM domains, or KEGG pathways in differentially expressed genes.
---

# AnoExpress — Enrichment analysis

Hypergeometric tests for over-representation of GO terms, PFAM domains, or KEGG pathways. All functions are importable from the top-level `anoexpress` namespace.

---

## Overview

All three enrichment functions share the same interface. You provide **either** a ranking function (`func`) to automatically select the top percentile of genes, **or** an explicit gene list (`gene_ids`) — not both.

The background set is all genes in the analysis that have at least one annotation of the relevant type.

Multiple testing correction uses Benjamini-Hochberg FDR (`statsmodels.stats.multitest.fdrcorrection`).

---

## `ae.go_hypergeometric()`

Test for over-representation of GO terms.

```python
ae.go_hypergeometric(
    analysis,            # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" | "fun"
    func=None,           # ranking function (e.g. np.nanmedian) — mutually exclusive with gene_ids
    gene_ids=None,       # explicit list of AGAP IDs — mutually exclusive with func
    percentile=0.05,     # fraction of top-ranked genes to use as the test set (when func is used)
    microarray=False,
    low_count_filter=None,
)
```

**Returns:** `pd.DataFrame` with columns `annotation` (GO term), `pval`, `padj`, `descriptions`.

**Example:**
```python
import numpy as np
import anoexpress as ae

result = ae.go_hypergeometric(
    analysis="gamb_colu_arab_fun",
    func=np.nanmedian,
    percentile=0.05,
)
result[result.padj < 0.05]
```

---

## `ae.pfam_hypergeometric()`

Test for over-representation of PFAM domains.

```python
ae.pfam_hypergeometric(
    analysis,
    func=None,
    gene_ids=None,
    percentile=0.05,
    microarray=False,
    low_count_filter=None,
)
```

**Returns:** `pd.DataFrame` with columns `annotation` (PFAM domain ID), `pval`, `padj`.

**Example:**
```python
result = ae.pfam_hypergeometric(
    analysis="gamb_colu_arab_fun",
    gene_ids=["AGAP004707", "AGAP002865", "AGAP009194"],
)
```

---

## `ae.kegg_hypergeometric()`

Test for over-representation of KEGG pathways.

```python
ae.kegg_hypergeometric(
    analysis,
    func=None,
    gene_ids=None,
    percentile=0.05,
    microarray=False,
    low_count_filter=None,
)
```

**Returns:** `pd.DataFrame` with columns `annotation` (KEGG pathway ID), `pval`, `padj`, `description`.

---

## Notes

- All three functions require a network connection; annotation files are downloaded from GitHub.
- PFAM and GO annotations are only available for *An. gambiae* (`AGAP` IDs). For `analysis="fun"`, use `gene_ids` containing AFUN orthologs where possible, or expect empty results.
- The `percentile` parameter is the fraction of the full ranked gene list to include as the "target" set (e.g. `0.05` = top 5%).
