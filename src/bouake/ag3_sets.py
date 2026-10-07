"""Which Ag3 sample sets may be used in the paper.

`ag3.sample_sets()["unrestricted_use"]` is False for sets still under their terms-of-use
embargo. Those sets can be read but must not enter a published analysis, so every Ag3 call in
this project is restricted to unrestricted sets (which also keeps each call small).
"""

import malariagen_data
import pandas as pd


def unrestricted_sample_sets(*, ag3: malariagen_data.Ag3, taxa: tuple[str, ...]) -> list[str]:
    """Sorted unrestricted sample sets that contain at least one sample of the given taxa."""
    sample_sets = ag3.sample_sets()
    allowed = set(sample_sets.loc[sample_sets["unrestricted_use"], "sample_set"])
    metadata = ag3.sample_metadata()
    with_taxa = set(metadata.loc[metadata["taxon"].isin(taxa), "sample_set"])
    return sorted(allowed & with_taxa)


def cap_per_cohort(*, metadata: "pd.DataFrame", cohort_column: str, max_per_cohort: int, seed: int = 42) -> "pd.DataFrame":
    """Random subsample of at most `max_per_cohort` samples per cohort; samples without a cohort are dropped.

    Adds a `cohort` column. Written without groupby.apply, which drops the grouping column in pandas 3.
    """
    metadata = metadata.assign(cohort=metadata[cohort_column].fillna("unassigned"))
    kept = metadata.sample(frac=1, random_state=seed).groupby("cohort", sort=False).head(max_per_cohort)
    return kept[kept["cohort"] != "unassigned"]
