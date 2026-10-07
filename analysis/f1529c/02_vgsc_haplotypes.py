"""A3 — haplotype background of F1529C in Ag3 An. coluzzii, by scaffold phasing.

The Ag3 phased haplotype set is biallelic by construction and holds only I1527T of the three
focal sites; V402L (multiallelic) and F1529C are absent. They are put back with the
scaffold-phasing method in selection-atlas (sweepclust/docs/SCAFFOLD-PHASING.md): each
site's unphased genotypes are phased onto the surrounding phased panel by a copying model.

    V402L[T]  2L:2391228 G>T   }  independent origins of V402L, phased separately
    V402L[C]  2L:2391228 G>C   }
    F1529C    2L:2429623 T>G
    I1527T    2L:2429617 T>C   already in the panel

A site is declined if more than 2% of its heterozygotes tie or it shows heterozygote excess
(F < -0.15), the refusals the method's validation supports.

Writes results/f1529c/vgsc_haplotypes.csv (one row per phased coluzzii haplotype) and
results/f1529c/vgsc_phasing_diagnostics.json.
"""

import json
import sys

import malariagen_data
import numpy as np
import pandas as pd

from bouake.ag3_sets import unrestricted_sample_sets
from bouake.paths import REPO

sys.path.insert(0, str(REPO.parent / "selection-atlas" / "sweepclust" / "src"))
from sweepatlas import data as sweep_data  # noqa: E402
from sweepatlas import scaffold  # noqa: E402

OUT = REPO / "results" / "f1529c"
# Wide enough that every target has a full 1,500-site scaffold window on both sides.
REGION = "2L:2380000-2440000"
SAMPLE_QUERY = "taxon == 'coluzzii'"
PANEL_SITES = {"I1527T": (2429617, "C")}
TARGETS = {"V402L_T": (2391228, "T"), "V402L_C": (2391228, "C"), "F1529C": (2429623, "G")}
MAX_TIE_FRACTION = 0.02
MIN_INBREEDING_F = -0.15


def inbreeding_f(*, dosage: np.ndarray) -> float:
    """1 - observed/expected heterozygosity; strongly negative flags a collapsed duplication."""
    called = dosage[dosage >= 0]
    p = called.mean() / 2.0
    expected = 2 * p * (1 - p)
    return float(1 - (called == 1).mean() / expected) if expected > 0 else float("nan")


def load_panel(*, ag3: malariagen_data.Ag3, sample_sets: list[str]):
    ds = ag3.haplotypes(region=REGION, analysis="gamb_colu", sample_sets=sample_sets,
                        sample_query=SAMPLE_QUERY)
    return (ds["call_genotype"].values, ds["variant_position"].values.astype(np.int64),
            ds["variant_allele"].values.astype(str), ds["sample_id"].values.astype(str))


def panel_haplotype_alleles(*, scaffold_calls, positions, alleles, position: int, alt: str) -> np.ndarray:
    """Derived-allele indicator per haplotype for a site already in the phased panel."""
    row = int(np.where(positions == position)[0][0])
    assert alleles[row].tolist() == ["T", "C"] if alt == "C" else alt in alleles[row]
    return (scaffold_calls[row] == 1).reshape(-1)


def phase_targets(*, ag3, scaffold_calls, positions, sample_ids, sample_sets) -> tuple[dict, dict]:
    """Phase each target onto the scaffold; returns haplotype vectors and diagnostics."""
    config = dict(sample_sets=sample_sets, sample_query=SAMPLE_QUERY)
    vectors, diagnostics = {}, {}
    for label, (position, alt) in TARGETS.items():
        site = sweep_data.fetch_unphased_site(ag3, "2L", position, alt, config, sample_ids)
        phased = scaffold.phase_site(scaffold_calls, positions, site.dosage, position)
        f = inbreeding_f(dosage=site.dosage)
        reasons = []
        if phased.tie_fraction > MAX_TIE_FRACTION:
            reasons.append(f"tie fraction {phased.tie_fraction:.3f}")
        if f < MIN_INBREEDING_F:
            reasons.append(f"heterozygote excess F = {f:.3f}")
        diagnostics[label] = dict(phased.summary(), ref=site.ref_allele, alt=alt,
                                  inbreeding_f=round(f, 4), declined_because=reasons)
        print(label, {k: diagnostics[label][k] for k in
                      ("n_het", "alt_freq", "tie_fraction", "converged", "inbreeding_f")},
              "DECLINED: " + "; ".join(reasons) if reasons else "")
        if not reasons:
            vectors[label] = phased.hap_allele
    return vectors, diagnostics


if __name__ == "__main__":
    ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
    metadata = ag3.sample_metadata()
    sample_sets = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii",))

    scaffold_calls, positions, alleles, sample_ids = load_panel(ag3=ag3, sample_sets=sample_sets)
    vectors, diagnostics = phase_targets(ag3=ag3, scaffold_calls=scaffold_calls, positions=positions,
                                         sample_ids=sample_ids, sample_sets=sample_sets)
    for label, (position, alt) in PANEL_SITES.items():
        vectors[label] = panel_haplotype_alleles(scaffold_calls=scaffold_calls, positions=positions,
                                                 alleles=alleles, position=position, alt=alt)

    table = pd.DataFrame({"sample_id": np.repeat(sample_ids, 2),
                          "haplotype": np.tile([0, 1], len(sample_ids))})
    for label, vector in vectors.items():
        table[label] = vector
    table = table.merge(metadata[["sample_id", "country", "admin1_iso", "year"]], on="sample_id")
    table.to_csv(OUT / "vgsc_haplotypes.csv", index=False)
    (OUT / "vgsc_phasing_diagnostics.json").write_text(json.dumps(diagnostics, indent=1, default=float))
    print(len(table), "haplotypes;", {k: int((table[k] == 1).sum()) for k in vectors})
