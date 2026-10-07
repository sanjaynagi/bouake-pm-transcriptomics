"""A1/A2 — Vgsc V402L, I1527T and F1529C allele frequencies across Ag3 cohorts.

The three focal codons are read directly from the Ag3 SNP calls (the library's
aa_allele_frequencies_advanced currently fails under the installed xarray/pandas). Each change
is identified by nucleotide, never by allele index:

    V402L   2L:2391228  G>T and G>C   two independent substitutions, both giving Leu
                                      (G>A gives V402I and is not counted)
    I1527T  2L:2429617  T>C
    F1529C  2L:2429623  T>G

Only unrestricted (not embargoed) sample sets that contain An. coluzzii or An. gambiae are loaded.

Writes
  results/f1529c/vgsc_sample_genotypes.csv   per sample: metadata + alt-allele copies at each site
  results/f1529c/vgsc_aa_frequencies.csv     per taxon x admin1 x year cohort (>= 10 samples)
"""


import malariagen_data
import numpy as np
import pandas as pd

from bouake.ag3_sets import unrestricted_sample_sets
from bouake.paths import REPO

OUT = REPO / "results" / "f1529c"
SITES = {
    "V402L": (2391228, "G", ("T", "C")),
    "I1527T": (2429617, "T", ("C",)),
    "F1529C": (2429623, "T", ("G",)),
}
METADATA_COLUMNS = ["sample_id", "taxon", "country", "admin1_iso", "admin1_name", "year", "sample_set"]
MIN_COHORT_SIZE = 10
TAXA = ("coluzzii", "gambiae")
SAMPLE_QUERY = "taxon in ['coluzzii', 'gambiae']"


def alt_copies(*, ag3: malariagen_data.Ag3, change: str, sample_sets: list[str]) -> pd.DataFrame:
    """Derived-allele copy number (0, 1, 2; -1 if uncalled) per sample at one focal site."""
    position, ref, alts = SITES[change]
    ds = ag3.snp_calls(region=f"2L:{position}-{position}", sample_sets=sample_sets,
                       sample_query=SAMPLE_QUERY, site_mask=None)
    alleles = [a.decode() if isinstance(a, bytes) else str(a) for a in ds["variant_allele"].values[0]]
    assert alleles[0] == ref and all(a in alleles for a in alts), (change, alleles)
    genotype = ds["call_genotype"].values[0]
    copies = np.isin(genotype, [alleles.index(a) for a in alts]).sum(axis=1).astype(int)
    copies[(genotype < 0).any(axis=1)] = -1
    print(change, "alleles", alleles, "derived", alts, flush=True)
    return pd.DataFrame({"sample_id": ds["sample_id"].values, change: copies})


def sample_genotypes(*, ag3: malariagen_data.Ag3) -> pd.DataFrame:
    """Metadata joined to derived-allele copy number at the three sites."""
    sample_sets = unrestricted_sample_sets(ag3=ag3, taxa=TAXA)
    metadata = ag3.sample_metadata()
    metadata = metadata[metadata["taxon"].isin(TAXA) & metadata["sample_set"].isin(sample_sets)]
    table = metadata[METADATA_COLUMNS]
    for change in SITES:
        table = table.merge(alt_copies(ag3=ag3, change=change, sample_sets=sample_sets),
                            on="sample_id", how="inner")
    return table


def cohort_frequencies(*, genotypes: pd.DataFrame) -> pd.DataFrame:
    """Per cohort allele frequency of each focal change (called samples only)."""
    rows = []
    for (taxon, country, admin1, year), cohort in genotypes.groupby(
        ["taxon", "country", "admin1_iso", "year"], dropna=False
    ):
        if len(cohort) < MIN_COHORT_SIZE:
            continue
        row = dict(taxon=taxon, country=country, admin1_iso=admin1, year=year, n_samples=len(cohort))
        for change in SITES:
            called = cohort.loc[cohort[change] >= 0, change]
            row[f"n_{change}"] = 2 * len(called)
            row[f"count_{change}"] = int(called.sum())
            row[f"frq_{change}"] = called.sum() / (2 * len(called)) if len(called) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
    genotypes = sample_genotypes(ag3=ag3)
    genotypes.to_csv(OUT / "vgsc_sample_genotypes.csv", index=False)
    frequencies = cohort_frequencies(genotypes=genotypes)
    frequencies.to_csv(OUT / "vgsc_aa_frequencies.csv", index=False)
    print(frequencies.shape)
