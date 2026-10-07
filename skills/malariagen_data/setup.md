---
name: malariagen-data-setup
description: API reference for setting up malariagen_data API objects (Ag3, Af1, Amin1, Adir1) and basic data access — releases, sample_sets, lookup_release, lookup_study. Use when initialising the API or discovering what data is available.
---

# malariagen_data — Setup and basic data access

## Initialising the API

All four mosquito APIs share the same constructor signature:

```python
import malariagen_data

ag3 = malariagen_data.Ag3()
af1 = malariagen_data.Af1()
amin1 = malariagen_data.Amin1()
adir1 = malariagen_data.Adir1()
```

### Constructor parameters (Ag3 shown; others are similar)

```python
malariagen_data.Ag3(
    url=None,                          # custom GCS URL (default: MalariaGEN public bucket)
    results_cache=None,                # path to local directory for caching computation results
    bokeh_output_notebook=True,        # call bokeh.output_notebook() on init
    log=sys.stdout,                    # log destination
    debug=False,                       # verbose debug output
    show_progress=None,                # show tqdm progress bars (None = auto-detect)
    check_location=True,               # warn if not running in GCS-adjacent compute
    cohorts_analysis=None,             # cohort analysis version string (default: latest)
    aim_analysis=None,                 # AIM analysis version (Ag3 only; default: latest)
    site_filters_analysis=None,        # site filters analysis version (default: latest)
    unrestricted_use_only=False,       # filter to only unrestricted-use samples
    surveillance_use_only=False,       # filter to only surveillance-use samples
    pre=False,                         # include pre-release data
    **storage_options,                 # passed to fsspec (e.g. token=...)
)
```

**Tip:** Set `results_cache="/tmp/malariagen_cache"` to persist expensive computations (PCA, GWSS, etc.) across sessions.

---

## `releases`

```python
ag3.releases
```

**Returns:** `tuple[str, ...]` — all available data release identifiers, e.g. `("3.0", "3.1", "3.2")`.

---

## `sample_sets()`

```python
ag3.sample_sets(
    release=None,   # str — filter to a specific release; None returns all
)
```

**Returns:** `pd.DataFrame` with columns `sample_set`, `study_id`, `study_url`, `terms_of_use_expiry_date`, `terms_of_use_url`, `release`.

**Example:**
```python
df = ag3.sample_sets(release="3.0")
df.head()
```

---

## `lookup_release()`

```python
ag3.lookup_release(sample_set="AG1000G-AO")
```

**Returns:** `str` — the release identifier the given sample set belongs to.

---

## `lookup_study()`

```python
ag3.lookup_study(sample_set="AG1000G-AO")
```

**Returns:** `str` — the study identifier associated with the sample set.

---

## Notes on data access

- Data is streamed from Google Cloud Storage. Run computations inside GCS (e.g. Google Colab) for best performance.
- Some data releases are subject to **terms of use** including embargoes. Check `terms_of_use_expiry_date` in `sample_sets()`.
- The `pre=True` flag unlocks access to pre-release data — use with caution.
- `unrestricted_use_only=True` or `surveillance_use_only=True` automatically filter all downstream calls to the appropriate data subset.
