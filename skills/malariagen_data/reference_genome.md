---
name: malariagen-data-reference-genome
description: API reference for malariagen_data reference genome functions — contigs, genome_sequence(), genome_features(), plot_transcript(), plot_genes(). Use when accessing genomic coordinates, gene annotations, or plotting gene/transcript structure.
---

# malariagen_data — Reference genome

---

## `contigs`

Property listing the chromosome/contig identifiers for the reference genome.

```python
ag3.contigs
# e.g. ('2R', '2L', '3R', '3L', 'X', 'Mt')
```

**Returns:** `tuple[str, ...]`

---

## `genome_sequence()`

Access the reference genome as a numpy array of nucleotide characters.

```python
ag3.genome_sequence(
    region,          # str — contig or region, e.g. "2L" or "2L:1000000-2000000"
    inline_array=True,
    chunks="native",
)
```

**Returns:** `numpy.ndarray` of `dtype="|S1"` (single-byte character per base).

**Example:**
```python
seq = ag3.genome_sequence("2L:1_000_000-1_100_000")
# Convert to string: seq.tobytes().decode()
```

---

## `genome_features()`

Access gene/feature annotations (from GFF3).

```python
ag3.genome_features(
    region=None,      # str or Region — filter to a genomic region; None = whole genome
    attributes=("ID", "Parent", "Name", "description"),  # GFF attributes to include
)
```

**Returns:** `pd.DataFrame` — one row per GFF feature. Columns include `contig`, `source`, `type`, `start`, `end`, `score`, `strand`, `phase`, plus one column per attribute.

**Common feature types:** `"gene"`, `"mRNA"`, `"exon"`, `"CDS"`, `"three_prime_UTR"`, `"five_prime_UTR"`

**Example:**
```python
# Get all genes on chromosome 2L
df_genes = ag3.genome_features(region="2L", attributes=("ID", "Name", "description"))
df_genes = df_genes.query("type == 'gene'")
```

---

## `plot_transcript()`

Plot the exon/intron structure of a transcript (gene model).

```python
ag3.plot_transcript(
    transcript,              # str — transcript ID, e.g. "AGAP004707-RA"
    sizing_mode="stretch_width",
    width=800,
    height=120,
    show=True,
    output_backend="svg",
)
```

**Returns:** Bokeh figure.

**Notes:**
- A single gene often has multiple transcripts. Use `genome_features()` to find transcript IDs (features with `type == "mRNA"`).
- This function is typically used to add a gene model track below a genomic plot.

---

## `plot_genes()`

Plot all genes in a genomic region as a track.

```python
ag3.plot_genes(
    region,                  # str or Region — e.g. "2L:28_000_000-29_000_000"
    sizing_mode="stretch_width",
    width=800,
    height=120,
    show=True,
    output_backend="svg",
    x_range=None,            # Bokeh Range1d — shared X axis for linking plots
)
```

**Returns:** Bokeh figure.

**Example — linked gene model track:**
```python
fig_h12 = ag3.plot_h12_gwss(contig="2L", ..., show=False)
fig_genes = ag3.plot_genes(region="2L", x_range=fig_h12.x_range, show=False)
bokeh.plotting.show(bokeh.layouts.column(fig_h12, fig_genes))
```
