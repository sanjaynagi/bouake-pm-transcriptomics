---
name: malariagen-data-sample-metadata
description: API reference for malariagen_data sample metadata functions — sample_metadata(), count_samples(), cohorts(), lookup_sample(), add_extra_metadata(), wgs_data_catalog(), and sample map/bar plot functions. Use when working with sample information, cohort definitions, or geographic visualisation.
---

# malariagen_data — Sample metadata

---

## `sample_metadata()`

The primary entry point for working with sample information.

```python
ag3.sample_metadata(
    sample_sets=None,            # str or list[str] — sample set IDs or release tag
    sample_query=None,           # str — pandas query string, e.g. "country == 'Ghana'"
    sample_query_options=None,   # dict — extra kwargs passed to DataFrame.query()
    sample_indices=None,         # list[int] — mutually exclusive with sample_query
)
```

**Returns:** `pd.DataFrame` — one row per sample. Key columns include:

| Column | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `country` | Country of collection |
| `admin1_iso` / `admin1_name` | Admin level 1 |
| `admin2_name` | Admin level 2 |
| `year` / `month` | Collection date |
| `latitude` / `longitude` | Collection coordinates |
| `taxon` | Inferred species/taxon |
| `sex_call` | Inferred sex (`"F"` / `"M"`) |
| `study_id` | Contributing study |

**Notes:**
- `sample_query` and `sample_indices` are mutually exclusive.
- When a query returns 0 results on a non-empty dataset, a `UserWarning` is raised with fuzzy-match suggestions for column names and values.

**Example:**
```python
df = ag3.sample_metadata(
    sample_sets="3.0",
    sample_query="country == 'Ghana' and taxon == 'gambiae'",
)
```

---

## `count_samples()`

Pivot table of sample counts by space, time, and taxon.

```python
ag3.count_samples(
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    index=("country", "admin1_iso", "admin1_name", "admin2_name", "year"),
    columns="taxon",
)
```

**Returns:** `pd.DataFrame` — pivot table of sample counts.

---

## `cohorts()`

Return a cohort table, grouping samples by taxon, area (admin level), and time period.

```python
ag3.cohorts(
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    min_cohort_size=10,           # int — minimum samples to include a cohort
)
```

**Returns:** `pd.DataFrame` — one row per cohort, with columns including `cohort_id`, `cohort_label`, `cohort_size`, `taxon`, `country`, `admin1_iso`, `admin2_name`, `year`.

---

## `lookup_sample()`

Look up metadata for a single sample.

```python
ag3.lookup_sample(sample="VBS00001-4248a")
```

**Returns:** `pd.Series` with all metadata columns for the given sample.

---

## `add_extra_metadata()`

Attach additional user-supplied columns to all future `sample_metadata()` calls.

```python
ag3.add_extra_metadata(
    data,      # pd.DataFrame — must contain a join key column
    on="sample_id",   # str — column to join on
)
```

---

## `clear_extra_metadata()`

Remove all previously added extra metadata.

```python
ag3.clear_extra_metadata()
```

---

## `wgs_data_catalog()`

List all sequencing run accessions for a sample set.

```python
ag3.wgs_data_catalog(sample_set="AG1000G-AO")
```

**Returns:** `pd.DataFrame` with columns including `sample_id`, `run_accession`, `bam_url`.

---

## `plot_samples_bar()`

Bar chart of sample counts by a metadata column.

```python
ag3.plot_samples_bar(
    x,                          # str — metadata column for X axis (e.g. "country")
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    color=None,                 # str — metadata column to colour bars
    sort=True,
    title="Sample counts",
    # ... Bokeh plot params: sizing_mode, width, height, show, output_backend
)
```

**Returns:** Bokeh figure.

---

## `plot_samples_interactive_map()`

Interactive Leaflet map of sample locations.

```python
ag3.plot_samples_interactive_map(
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    min_cohort_size=10,
    # ... map display params
)
```

**Returns:** `ipyleaflet.Map`

---

## `plot_sample_location_mapbox()` / `plot_sample_location_geo()`

Plotly map visualisations of sample locations.

```python
ag3.plot_sample_location_mapbox(
    sample_sets=None,
    sample_query=None,
    sample_query_options=None,
    sample_indices=None,
    color="taxon",
    # ... Plotly params
)

ag3.plot_sample_location_geo(
    sample_sets=None,
    sample_query=None,
    # ... Plotly params
)
```

**Returns:** Plotly figure.

---

## `cross_metadata` (Ag3 only)

Property returning metadata for lab crosses included in the Ag3 data.

```python
df_crosses = ag3.cross_metadata
```

**Returns:** `pd.DataFrame`
