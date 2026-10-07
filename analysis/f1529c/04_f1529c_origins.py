"""A3 — one origin of F1529C or two?

F1529C sits on 1527T haplotypes carrying either of the two independent 402L substitutions
(script 03). If 1529C arose once, all carriers share an extended haplotype around codon 1529
whatever their 402L origin (402L is 38 kb away, so recombination can separate it from 1529);
if it arose separately on each background, carriers of different 402L origin will be no more
alike than non-carriers are.

Over windows of phased SNPs around codon 1529 (the 1527 and 1529 sites themselves excluded) the
pairwise Hamming distance is compared between:
    carrier - carrier              all 1529C haplotypes
    carrier[C] - carrier[T]        carriers of the two 402L origins
    carrier[C] - carrier[C], carrier[T] - carrier[T]   carriers of the same 402L origin
    control[C] - control[T]        non-carriers of the two 402L origins: whether the two 402L
                                   backgrounds are interchangeable across the whole 1527T lineage
    carrier - control              carriers against a random sample of non-carrier 1527T haplotypes
    control - control              non-carriers against each other
and the share of carrier pairs that are identical is reported. The haplotype x SNP allele matrix over the widest window, for every carrier plus 150
non-carriers per 402L background, is saved for the haplotype plot.

Haplotype-level bootstrap intervals (500 resamples) are added for the widest window.

Writes results/f1529c/f1529c_origins_summary.csv, f1529c_origins_haplotypes.npz and
f1529c_lineage_marker.csv.
"""


import malariagen_data
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist

from bouake.ag3_sets import unrestricted_sample_sets
from bouake.paths import REPO

OUT = REPO / "results" / "f1529c"
FOCAL_POSITION = 2429623
EXCLUDED_POSITIONS = (2429617, 2429623)
REGION = "2L:2400000-2460000"
WINDOWS_BP = (2_500, 5_000, 10_000, 25_000)
N_CONTROLS = 900
N_BOOTSTRAP = 500
N_TREE_CONTROLS = 150
RNG = np.random.default_rng(42)


def load_panel(*, ag3: malariagen_data.Ag3, sample_sets: list[str]):
    ds = ag3.haplotypes(region=REGION, analysis="gamb_colu", sample_sets=sample_sets,
                        sample_query="taxon == 'coluzzii'")
    return (ds["call_genotype"].values, ds["variant_position"].values.astype(np.int64),
            ds["sample_id"].values.astype(str))


def haplotype_matrix(*, calls: np.ndarray, sample_ids: np.ndarray, table: pd.DataFrame) -> np.ndarray:
    """Haplotype x site allele matrix (all panel sites) for the resolved 1527T haplotypes."""
    row_of = {sample: i for i, sample in enumerate(sample_ids)}
    return np.stack([calls[:, row_of[s], h] for s, h in
                     zip(table["sample_id"], table["haplotype"])]).astype(np.int8)


def pair_summary(*, distance: np.ndarray, same_set: bool, n_boot: int = 0) -> dict:
    """Median distance and share of identical pairs; optional haplotype-level bootstrap interval.

    The bootstrap resamples haplotypes (not pairs), because pairs share haplotypes and are not
    independent. In a within-set comparison a resampled haplotype is never paired with itself.
    """
    values = distance[np.triu_indices_from(distance, k=1)] if same_set else distance.ravel()
    result = dict(median=float(np.median(values)), share_identical=float((values == 0).mean()))
    if n_boot:
        identical = distance == 0
        shares = []
        for _ in range(n_boot):
            rows = RNG.integers(0, distance.shape[0], distance.shape[0])
            cols = RNG.integers(0, distance.shape[1], distance.shape[1])
            block = identical[np.ix_(rows, cols)]
            if same_set:
                block = block[rows[:, None] != cols[None, :]]
            shares.append(block.mean())
        result["ci_low"], result["ci_high"] = (float(x) for x in np.percentile(shares, [2.5, 97.5]))
    return result


def lineage_marker(*, positions, controls, matrix, origin, carriers_c, carriers_t) -> None:
    """The SNP that best separates non-carrier 402L[G>T] from 402L[G>C] haplotypes, and who carries it.

    The site is chosen using non-carriers only, so carriers are an independent test. If 1529C had
    arisen separately on the G>T lineage, G>T carriers would carry that lineage's allele; if it
    arose once on the G>C lineage they would not.
    """
    controls_c, controls_t = controls[origin[controls] == "C"], controls[origin[controls] == "T"]
    difference = np.abs(matrix[controls_t].mean(axis=0) - matrix[controls_c].mean(axis=0))
    site = int(np.argmax(difference))
    lineage_allele = int(matrix[controls_t, site].mean() > 0.5)       # the allele common in G>T controls
    rows = []
    for name, group in (("non-carriers, 402L G>T", controls_t), ("non-carriers, 402L G>C", controls_c),
                        ("F1529C carriers, 402L G>T", carriers_t), ("F1529C carriers, 402L G>C", carriers_c)):
        rows.append(dict(group=name, n_haplotypes=len(group), share_with_lineage_allele=float(
            (matrix[group, site] == lineage_allele).mean())))
    table = pd.DataFrame(rows).assign(position=int(positions[site]), offset_bp=int(positions[site] - FOCAL_POSITION))
    table.to_csv(OUT / "f1529c_lineage_marker.csv", index=False)
    print(table.to_string(index=False), flush=True)


if __name__ == "__main__":
    ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
    sample_sets = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii",))
    table = pd.read_csv(OUT / "f1529c_haplotype_background.csv").reset_index(drop=True)
    calls, positions, sample_ids = load_panel(ag3=ag3, sample_sets=sample_sets)
    matrix = haplotype_matrix(calls=calls, sample_ids=sample_ids, table=table)

    is_carrier = (table["f1529c"] == 1).to_numpy()
    origin = table["v402l_origin"].to_numpy()
    carriers_c = np.where(is_carrier & (origin == "C"))[0]
    carriers_t = np.where(is_carrier & (origin == "T"))[0]
    carriers = np.where(is_carrier)[0]
    controls = RNG.choice(np.where(~is_carrier)[0], size=N_CONTROLS, replace=False)
    controls_c, controls_t = controls[origin[controls] == "C"], controls[origin[controls] == "T"]

    rows = []
    for window in WINDOWS_BP:
        keep = np.where((np.abs(positions - FOCAL_POSITION) <= window)
                        & ~np.isin(positions, EXCLUDED_POSITIONS))[0]
        sub = matrix[:, keep]
        distance = lambda a, b: cdist(sub[a], sub[b], metric="hamming") * len(keep)
        n_boot = N_BOOTSTRAP if window == WINDOWS_BP[-1] else 0
        comparisons = {
            "carrier - carrier": pair_summary(distance=distance(carriers, carriers), same_set=True, n_boot=n_boot),
            "carrier[C] - carrier[T]": pair_summary(distance=distance(carriers_c, carriers_t), same_set=False, n_boot=n_boot),
            "carrier[C] - carrier[C]": pair_summary(distance=distance(carriers_c, carriers_c), same_set=True, n_boot=n_boot),
            "carrier[T] - carrier[T]": pair_summary(distance=distance(carriers_t, carriers_t), same_set=True, n_boot=n_boot),
            "control[C] - control[T]": pair_summary(distance=distance(controls_c, controls_t), same_set=False, n_boot=n_boot),
            "carrier - control": pair_summary(distance=distance(carriers, controls), same_set=False, n_boot=n_boot),
            "control - control": pair_summary(distance=distance(controls, controls), same_set=True, n_boot=n_boot),
        }
        for name, result in comparisons.items():
            rows.append(dict(window_bp=window, n_sites=len(keep), comparison=name, **result))
        print(pd.DataFrame([r for r in rows if r["window_bp"] == window]).to_string(index=False), flush=True)
        if window == WINDOWS_BP[-1]:
            # Subset for the origins tree: every carrier plus equal numbers of non-carriers on
            # each 402L background, so neither background is drowned out.
            others = np.where(~is_carrier)[0]
            tree_controls = np.concatenate([
                RNG.choice(others[origin[others] == background],
                           size=min(N_TREE_CONTROLS, int((origin[others] == background).sum())), replace=False)
                for background in ("C", "T")])
            subset = np.sort(np.concatenate([carriers, tree_controls]))
            np.savez_compressed(OUT / "f1529c_origins_haplotypes.npz", subset=subset,
                                alleles=sub[subset], positions=positions[keep])
            lineage_marker(positions=positions[keep], controls=tree_controls,
                           matrix=sub, origin=origin, carriers_c=carriers_c, carriers_t=carriers_t)
    pd.DataFrame(rows).to_csv(OUT / "f1529c_origins_summary.csv", index=False)
