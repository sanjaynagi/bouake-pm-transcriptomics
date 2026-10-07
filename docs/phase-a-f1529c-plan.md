
# Phase A — Vgsc F1529C triple-mutant plan

## Problem

The paper reports a novel mutation, *Vgsc*-F1529C, in *An. coluzzii* from Bouaké (Côte d'Ivoire), equivalent to the well-characterised *Ae. aegypti* F1534C pyrethroid-resistance mutation (Musca codon numbering). Current manuscript evidence rests on:

1. A single IGV snapshot of one sample (`coluzziiPM3`) suggesting 1529C is linked to I1527T on the same reads.
2. Mean allele frequencies from RNA-Seq read balance.

This is insufficient for a claim of a triple mutant (402L/1527T/1529C) or for conclusions about geographic origin. Two things are unknown:

- Is F1529C present in wild *An. coluzzii* populations elsewhere, and if so, at what frequencies?
- Does F1529C sit on the 402L/1527T haplotype consistently at population scale, or only in the single Bouaké read pileup?

## Proposed solution

### A1. Geographic frequency survey (aa-level)

Use `ag3.aa_allele_frequencies_advanced(transcript="AGAP004707-RD", area_by="admin1_iso", period_by="year", sample_query="taxon == 'coluzzii'")` across all Ag3 releases, with `min_cohort_size=10`. Output is an xarray `Dataset` covering every coluzzii admin1-year cohort.

Extract F1529C, I1527T, and V402L rows. Report per-cohort frequencies and sample sizes. Compare against the separately computed *An. gambiae* distribution (expected to lack 1529C). Expected result: F1529C either absent elsewhere (→ Bouaké origin claim supportable) or present in a clustered subset of West African coluzzii cohorts (→ regional spread).

Complement with `ag3.plot_frequencies_interactive_map(ds, variant="1529C", ...)` for a geographic visual.

### A2. Co-occurrence analysis (aa-level)

Using the same `Dataset`, compute Pearson correlation of per-cohort frequencies between the three focal amino-acid changes. If 1529C arose on the 402L/1527T haplotype, its cohort-level frequency should be bounded by and correlated with both. A scatter (frq_1529C vs frq_1527T, coloured by frq_402L) makes the population-scale linkage visible.

### A3. Haplotype-level evidence

Pull phased haplotype data for the *Vgsc* region:

```python
ds_haps = ag3.haplotypes(
    region="2L:2_390_000-2_440_000",
    analysis="gamb_colu",
    sample_query="taxon == 'coluzzii'",
)
```

For each haplotype, extract the allele state at positions 402, 1527, 1529 (from `haplotype_sites()` + gene coordinates). Tabulate the joint haplotype (e.g. 402V/1527I/1529F wildtype vs 402L/1527T/1529F vs 402L/1527T/1529C) and its frequency across cohorts.

This replaces read-level IGV evidence with phased-haplotype evidence across hundreds of wild genomes. Expected: 1529C occurs exclusively on the 402L/1527T background.

If the Bouaké samples themselves lack phased haplotypes in Ag3 (RNA-Seq / low-coverage), rely on the wild-population phased data to establish the haplotype identity, and cross-reference the Bouaké VOI allele balance table (`results/variantsOfInterest/csvs/`) for within-sample co-occurrence.

### A4–A6. Outputs

- **Figure 4 panel (new)** — admin1-level choropleth or point map of F1529C frequency in coluzzii, with 402L/1527T baseline as context.
- **Figure 4 panel or Supplementary** — haplotype frequency stacked bar plot split by 402/1527/1529 status.
- **Supplementary table** — per-cohort frequencies for all three positions, plus sample counts.

## Risks and tradeoffs

- **Phasing coverage** — Ag3 phased haplotypes may not cover every coluzzii cohort; some cohorts will drop out. Document the set used.
- **Position numbering** — Ag3 uses AgamP4 coordinates; triple-check position 1529 corresponds to the correct codon in transcript AGAP004707-RD (the standard *Vgsc* transcript). Cross-reference against Clarkson et al. 2021 supplementary for 402L/1527T numbering.
- **Sample set freshness** — Bouaké 2019 samples may not be in Ag3. The haplotype analysis is therefore about establishing what the 402L/1527T/1529C haplotype *looks like* in the wild Ag3 catalogue; linking Bouaké specifically relies on the existing RNA-Seq allele-balance data already in the paper.
- **Null-result risk** — F1529C may be absent from Ag3. That would still be a positive finding (novel to Bouaké / CdI) and a clean geographic panel.

## Success criteria

- [ ] Per-cohort frequency table for V402L, I1527T, F1529C covering all coluzzii Ag3 cohorts ≥10 samples.
- [ ] Either (a) clear evidence F1529C is not present elsewhere in Ag3 coluzzii, or (b) a defined geographic distribution.
- [ ] Phased-haplotype evidence that 1529C sits on the 402L/1527T background (or explicit null if 1529C absent from Ag3).
- [ ] One publication-quality figure panel (geographic) + one supplementary figure (haplotype) + one supplementary table.
- [ ] Manuscript paragraph draft under "Variants of interest" replacing the IGV-only evidence with population-scale evidence.

## Deliverables

- `notebooks/f1529c/01_vgsc_aa_frequencies.ipynb`
- `notebooks/f1529c/02_vgsc_haplotypes.ipynb`
- `notebooks/f1529c/03_figure_panel.ipynb`
- `results/f1529c/vgsc_aa_frequencies.csv`
- `results/f1529c/vgsc_haplotype_table.csv`
- `figures_ms/f1529c_geo_map.{svg,png}`
- `figures_ms/f1529c_haplotype.{svg,png}`
