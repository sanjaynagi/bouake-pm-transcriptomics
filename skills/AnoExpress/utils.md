---
name: anoexpress-utils
description: API reference for AnoExpress utility functions — ae.resolve_gene_id(), ae.load_gff(). Use when resolving gene identifiers from genomic coordinates or file paths, or when loading GFF genome feature data.
---

# AnoExpress — Utilities

Helper functions used internally and available in the public API. Importable from the top-level `anoexpress` namespace.

---

## `ae.resolve_gene_id()`

Normalise a gene identifier into a list of AGAP/AFUN IDs. Handles four input formats:

| Input | Action |
|---|---|
| Genomic span `"2L:500000-600000"` | Queries the GFF for overlapping protein-coding genes |
| `.tsv` / `.txt` file path | Reads first column as gene IDs |
| `.csv` file path | Reads first column as gene IDs |
| `.xlsx` file path | Reads first column as gene IDs |
| Already a list | Returned as-is |

```python
ae.resolve_gene_id(
    gene_id,              # str or list
    analysis,             # used to validate genomic-span queries (not supported for "fun")
    gff_method='vectorbase',  # "vectorbase" | "malariagen_data"
    use_cache=True,
)
```

**Returns:** `list` of gene ID strings.

**Example:**
```python
# genes overlapping a genomic window
genes = ae.resolve_gene_id("2L:28460000-28560000", analysis="gamb_colu")

# from a file
genes = ae.resolve_gene_id("/path/to/my_genes.tsv", analysis="gamb_colu")
```

---

## `ae.load_gff()`

Load the *Anopheles gambiae* genome feature file (GFF) as a DataFrame. Data is cached at `~/.cache/anoexpress/` to avoid repeated downloads.

```python
ae.load_gff(
    method='vectorbase',   # "vectorbase" | "malariagen_data"
    override_type=None,    # override the default feature type filter
    query=None,            # pandas query string to filter rows
    use_cache=True,
)
```

**Returns:** `pd.DataFrame` with columns `contig`, `source`, `type`, `start`, `end`, `strand`, `attributes`, `GeneID` (and others from the raw GFF).

**Notes:**
- `method="vectorbase"`: downloads from VectorBase release 68, filters to `protein_coding_gene` by default. Chromosome arm coordinates are adjusted so that 2L/3L positions are offset onto combined `2RL`/`3RL` contigs (matching MalariaGEN conventions).
- `method="malariagen_data"`: uses the `malariagen_data` package to load genome features; filters to `gene` type by default.
- Cache files are stored as Parquet under `~/.cache/anoexpress/`.

**Example:**
```python
gff = ae.load_gff(
    method="vectorbase",
    query="contig == '2RL' and start > 28000000 and end < 29000000",
)
```

---

## Internal helper

### `_gene_ids_from_annotation(gene_annot_df, annotation)`

Extract gene IDs matching a GO term or PFAM domain from the combined annotation DataFrame returned by `ae.load_annotations()`. Used internally by `plot_gene_family_expression()` and `load_candidates()`.

- Strings starting with `"GO:"` are matched against the `GO_terms` column.
- All other strings are matched against the `domain` (PFAM) column.
- Accepts a single string or a list.
