"""A3 — what haplotype background does F1529C sit on?

F1529C could not be phased onto the scaffold (script 02 declines it: 29% of its heterozygotes
tie), so it is placed by logic instead. Every F1529C carrier also carries I1527T (script 01), so
1529C can only be on a 1527T haplotype. In a sample with n phased 1527T haplotypes:

    1529C copies == 0   no haplotype carries 1529C
    1529C copies == n   every 1527T haplotype carries 1529C
    otherwise           ambiguous, excluded

402L status of each haplotype is read from the phased V402L[T] and V402L[C] alleles (script 02),
so the two independent 402L origins are both counted. Writes
results/f1529c/f1529c_haplotype_background.csv and prints the summary used in the paper.
"""


import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from bouake.paths import RESULTS

OUT = RESULTS / "f1529c"


def f1529c_state(*, copies_1529c: int, copies_1527t: int) -> float:
    """1.0 / 0.0 for every 1527T haplotype of a sample, or NaN when the sample is ambiguous."""
    if copies_1529c == 0:
        return 0.0
    if copies_1529c == copies_1527t:
        return 1.0
    return np.nan


def background_table() -> pd.DataFrame:
    haplotypes = pd.read_csv(OUT / "vgsc_haplotypes.csv")
    genotypes = pd.read_csv(OUT / "vgsc_sample_genotypes.csv")[["sample_id", "F1529C"]]
    haplotypes = haplotypes.merge(genotypes, on="sample_id", how="left")
    n_1527t = haplotypes.groupby("sample_id")["I1527T"].transform("sum")
    haplotypes = haplotypes[haplotypes["I1527T"] == 1].copy()
    haplotypes["n_1527T_in_sample"] = n_1527t[haplotypes.index]
    haplotypes["f1529c"] = [
        f1529c_state(copies_1529c=int(c), copies_1527t=int(n)) if c >= 0 else np.nan
        for c, n in zip(haplotypes["F1529C"], haplotypes["n_1527T_in_sample"])
    ]
    haplotypes["v402l_origin"] = np.select(
        [haplotypes["V402L_T"] == 1, haplotypes["V402L_C"] == 1], ["T", "C"], default="none")
    return haplotypes.dropna(subset=["f1529c"])


if __name__ == "__main__":
    table = background_table()
    table.to_csv(OUT / "f1529c_haplotype_background.csv", index=False)
    carrying, other = table[table["f1529c"] == 1], table[table["f1529c"] == 0]
    print(f"1527T haplotypes resolved: {len(table)}; with 1529C: {len(carrying)}")
    print("\n402L origin on 1527T haplotypes, with vs without 1529C")
    print(pd.crosstab(table["f1529c"], table["v402l_origin"], margins=True))
    print(pd.crosstab(table["f1529c"], table["v402l_origin"], normalize="index").round(3))
    with_402l = lambda frame: int((frame["v402l_origin"] != "none").sum())
    contingency = [[with_402l(carrying), len(carrying) - with_402l(carrying)],
                   [with_402l(other), len(other) - with_402l(other)]]
    print("\nFisher exact (402L present vs absent, 1529C vs not):", contingency,
          fisher_exact(contingency))
    print("\n1529C haplotypes by country and 402L origin")
    print(pd.crosstab(carrying["country"], carrying["v402l_origin"], margins=True))
