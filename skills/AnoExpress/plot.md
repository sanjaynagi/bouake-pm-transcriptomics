---
name: anoexpress-plot
description: API reference for AnoExpress plotting functions — ae.plot_gene_expression(), ae.plot_gene_family_expression(), ae.plot_heatmap(), ae.plot_contig_expression_track(), ae.plot_contig_expression(). Use when visualising gene expression from the AnoExpress meta-analysis.
---

# AnoExpress — Plotting

All plot functions return interactive figures (Plotly or Bokeh). They are importable from the top-level `anoexpress` namespace.

---

## `ae.plot_gene_expression()`

Plot fold changes (and raw counts) for one or more genes across all experiments. Produces a two-panel interactive Plotly figure: fold-change strip/box plot on the left, raw counts on the right.

```python
ae.plot_gene_expression(
    gene_id,                      # str | list | file path
    analysis="gamb_colu_arab_fun",
    microarray=False,
    sample_query=None,            # pandas query to subset comparisons
    title=None,
    plot_type='strip',            # "strip" | "boxplot"
    sort_by='agap',               # "median" | "mean" | "agap" | "position" | None
    gff_method='vectorbase',
    pvalue_filter=None,           # float — remove non-significant FCs from plot
    width=1600,
    height=None,                  # auto-sized if None
    save_html=None,               # path to save as .html
)
```

**Returns:** `plotly.graph_objects.Figure`

**Example:**
```python
fig = ae.plot_gene_expression(
    gene_id=["AGAP004707", "AGAP002865"],
    analysis="gamb_colu_arab_fun",
    pvalue_filter=0.05,
    sort_by="median",
)
fig.show()
```

---

## `ae.plot_gene_family_expression()`

Like `plot_gene_expression()` but for a whole gene family, looked up by GO term or PFAM domain.

```python
ae.plot_gene_family_expression(
    gene_identifier,              # GO term (e.g. "GO:0004364") or PFAM domain (e.g. "GST_N")
    analysis,                     # "gamb_colu" | "gamb_colu_arab" | "gamb_colu_arab_fun" (not "fun")
    title,
    microarray=False,
    plot_type='strip',            # "strip" | "boxplot"
    sort_by='median',
    gff_method='vectorbase',
    width=1600,
    height=None,
)
```

**Returns:** `plotly.graph_objects.Figure`

**Note:** `analysis="fun"` is not supported because GO/PFAM annotations are only available for *An. gambiae*.

**Example:**
```python
fig = ae.plot_gene_family_expression(
    gene_identifier="GST_N",
    analysis="gamb_colu_arab_fun",
    title="Glutathione S-transferases",
)
fig.show()
```

---

## `ae.plot_heatmap()`

Seaborn clustermap of fold changes for selected genes.

```python
ae.plot_heatmap(
    analysis,
    gene_id=None,            # explicit gene list; if None, uses load_candidates()
    query_annotation=None,   # GO term or PFAM domain filter
    query_func=np.nanmedian, # function to rank genes when gene_id is None
    query_fc=None,           # minimum fold change threshold
    query_name='median',
    cmap=None,               # seaborn colormap name
    cbar_pos=None,           # "left" | "right" | "top" | "bottom" | None
    figsize=None,            # (width, height) tuple
)
```

**Returns:** `None` (renders inline). The underlying `seaborn.ClusterGrid` is not returned.

**Notes:**
- If `gene_id` is provided, `query_annotation` and `query_fc` are ignored.
- If neither `gene_id` nor `query_annotation` is given, the top genes ranked by `query_func` are plotted.
- Columns are individual comparisons; rows are genes labelled `AGAP | GeneName`.

**Example:**
```python
import numpy as np
ae.plot_heatmap(
    analysis="gamb_colu_arab_fun",
    query_func=np.nanmedian,
    query_fc=2,
    query_name="median",
)
```

---

## `ae.plot_contig_expression_track()`

Genome-wide expression scan for a single chromosome arm. Plots raw fold changes per gene (circles) and a sliding-window median (line) using Bokeh.

```python
ae.plot_contig_expression_track(
    contig,              # "2L" | "2R" | "3L" | "3R" | "X" | "2RL" | "3RL"
    analysis="gamb_colu_arab_fun",
    data_type='fcs',     # "fcs" | "log2counts"
    microarray=False,
    pvalue_filter=None,
    size=10,             # window size (genes) for moving average
    step=5,              # step size (genes) for moving average
    title=None,
    width=800,
    height=600,
    palette=None,        # list of 4 colours for [coluzzii, gambiae, arabiensis, funestus]
    sizing_mode='stretch_width',
    x_range=None,        # bokeh Range1d
    y_range=None,
    show=False,
)
```

**Returns:** `bokeh.plotting.Figure`

---

## `ae.plot_contig_expression()`

Like `plot_contig_expression_track()` but adds a gene-model track beneath using `malariagen_data`. Requires `malariagen_data` to be installed and configured.

```python
ae.plot_contig_expression(
    contig,
    analysis,
    data_type='fcs',
    microarray=False,
    size=10,
    step=5,
    pvalue_filter=None,
    palette=None,
    y_range=(-10, 15),
    height=400,
    width=600,
    title=None,
    show=False,
)
```

**Returns:** `bokeh.layouts.GridPlot`
