---
name: malariagen-data-aim
description: API reference for malariagen_data ancestry-informative markers (AIMs) — aim_ids, aim_variants(), aim_calls(), plot_aim_heatmap(). Ag3 only. Use when classifying samples by species using AIM genotypes or visualising admixture.
---

# malariagen_data — AIM data (Ag3 only)

Ancestry-informative markers (AIMs) are SNPs that are highly differentiated between species and are used to classify samples into *An. gambiae*, *An. coluzzii*, and *An. arabiensis*.

---

## `aim_ids`

Property listing available AIM set identifiers.

```python
ag3.aim_ids
# e.g. ('gambcolu_vs_arab', 'gamb_vs_colu')
```

---

## `aim_variants()`

Access the variant sites for a given AIM set.

```python
ag3.aim_variants(aims="gambcolu_vs_arab")
```

**Returns:** `xarray.Dataset` with variant position and allele information.

---

## `aim_calls()`

Access AIM genotype calls for a set of samples.

```python
ag3.aim_calls(
    aims,                        # str — AIM set ID, e.g. "gambcolu_vs_arab"
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    inline_array=True,
    chunks="native",
)
```

**Returns:** `xarray.Dataset`

---

## `plot_aim_heatmap()`

Heatmap of AIM genotypes across samples, useful for visualising admixture.

```python
ag3.plot_aim_heatmap(
    aims,
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    sort_by_aims=True,
    sizing_mode="stretch_width",
    width=800,
    height=600,
    show=True,
    output_backend="svg",
)
```

**Returns:** Bokeh figure.
