---
name: anoexpress
description: Overview and index for the AnoExpress Python package — Anopheles gene expression in insecticide resistance studies. Use this skill to understand what the package does and navigate to the right API reference.
---

# AnoExpress

**Anopheles gene expression in resistance studies**

A Python package and meta-analysis of RNA-Seq (and microarray via IR-Tex) studies investigating insecticide resistance in *Anopheles gambiae s.l* and *Anopheles funestus*. Gene expression data is accessed remotely from GitHub; no local data files are required.

## Installation

```bash
pip install anoexpress
```

```python
import anoexpress as ae
```

## Key concepts

### `analysis` parameter

All major functions require an `analysis` argument that selects the set of species included. Analyses with more species contain fewer genes (due to ortholog finding):

| Value | Species included |
|---|---|
| `"gamb_colu"` | *An. gambiae*, *An. coluzzii* |
| `"gamb_colu_arab"` | above + *An. arabiensis* |
| `"gamb_colu_arab_fun"` | above + *An. funestus* |
| `"fun"` | *An. funestus* only |

### `data_type` parameter

| Value | Description |
|---|---|
| `"fcs"` | Log2 fold changes (resistant vs susceptible) |
| `"pvals"` | Adjusted p-values from differential expression |
| `"log2counts"` | Normalised log2 read counts |

### Gene IDs

- *An. gambiae* genes use `AGAP` identifiers (e.g. `"AGAP004707"`)
- *An. funestus* genes use `AFUN` identifiers
- Genomic spans are accepted as `"2L:500000-600000"`
- A path to a `.tsv`, `.txt`, `.csv`, or `.xlsx` file (gene IDs in first column) is also accepted

### `microarray` parameter

When `True`, IR-Tex microarray data is merged with the RNA-Seq data for a broader comparison set.

---

## Skills index

| Skill file | Contents |
|---|---|
| [data.md](data.md) | Loading expression data and metadata |
| [plot.md](plot.md) | Visualisation functions |
| [candidates.md](candidates.md) | Ranking and filtering candidate genes |
| [enrichment.md](enrichment.md) | GO, PFAM and KEGG hypergeometric enrichment |
| [utils.md](utils.md) | Utility / helper functions |
