# Analysis code

Shared code lives in the `bouake` package (`src/bouake/`, on `PYTHONPATH` through pixi):
`bouake.paths` (project locations) and `bouake.ag3_sets` (restricts every Ag3 call to sample sets
without a use restriction, `unrestricted_use`; never load Ag3 without it).

Run everything with `pixi run python <script>`; scripts write to `results/<topic>/` and the
manuscript figure scripts in `figures_ms/` read from there.

## abc_obp/ — ABC transporters, odorant-binding proteins and the ABCH1 locus
| Script | Produces |
|---|---|
| `setup.py` | shared gene families (strict ABC = PF00005 or `Abc*`), DE tables, sweep overlap; caches in `results/abc_obp/cache/` |
| `b9_abc_tables.py` | `bouake_de_abc_*.csv`, `sweep_overlap_abc_gambiae.csv` |
| `a_sweeps.py`, `a2_collapse.py`, `b_stats.py`, `c_obp_candidate.py` | sweep-overlap, clustering and OBP statistics |
| `h12_abch1.py`, `extra_stats.py`, `pbs.py`, `cnv_h1x.py` | H12, iHS/G123, PBS, copy number and H1X at ABCH1 |
| `hapclust_both.py` | ABCH1 haplotype clustering and NJ tree (both species) used in Figure 4 |
| `hapclust_abch1.py`, `temporal_tags.py` | coluzzii-only clustering and tag-SNP time course (exploratory, not in the paper) |
| `plot_clusters.py`, `plot_pbs.py` | exploratory plots |
| `legacy_notebooks/` | early notebooks that still use the old, broader ABC definition; superseded by the scripts above |

## f1529c/ — Vgsc F1529C
`01`–`04` run in order (frequency survey, phased haplotypes, haplotype background, single-origin test).
`03`/`04` feed Supplementary Figure 2; the survey (`01`) and cohort figures are kept for a separate Ag3 paper.

## Manuscript figures (`figures_ms/`)
`figure_gsea.py` (Fig. 2), `figure_variants.py` (Fig. 3), `figure_abch1.py` (Fig. 4), `figure_supp_f1529c.py` (Supp. Fig. 2),
with `figure_f1529c*.py` for the future Ag3 paper.

## rnaseq/ — transcriptomic results (reads processed by RNA-Seq-Pop, see `workflow/`)
| Paper result | Code | Input in this repository |
|---|---|---|
| DE values quoted in the text | `rnaseq/01_de_numbers.py` | `results/BouakePM_diffexp.xlsx` |
| Fig. 2 enrichment | `figures_ms/figure_gsea.py` | `figures_ms/gsea-bouake.xlsx` |
| Fig. 3 variants of interest | `figures_ms/figure_variants.py` | `results/variantsOfInterest/csvs/` |
| ABC/OBP, sweeps, ABCH1 (Fig. 4) | `abc_obp/`, `figures_ms/figure_abch1.py` | Ag3 via `malariagen_data` |
| F1529C (Supp. Fig. 2) | `f1529c/`, `figures_ms/figure_supp_f1529c.py` | Ag3 and Fig. 3 read counts |

## Not yet reproducible from this repository
These results come from RNA-Seq-Pop outputs or files that are not in the repo; the old exploratory notebooks
(`figures_ms/*.ipynb`, `results/notebooks/`) point at paths on other machines.
- Read and alignment totals (1.463 billion reads; 85.71% alignment) — Kallisto/fastp summaries.
- PCA and volcano plots (Fig. 1B, C) — normalised counts and gene-level DE tables (`old_results/genediff`).
- Bioassay mortality (Fig. 1A) — `resources/bioassay_data.xlsx`.
- Ancestry (AIM) proportions — per-sample AIM table; the values in the manuscript are typed into `figure_variants.py`.
- Karyotype frequencies and Welch tests — per-sample karyotype tables. The rendered notebook
  `results/notebooks/karyotype.ipynb` lists 21 samples (three coluzzii control and survivor samples) and does
  not reproduce the P values in the paper (2Rb coluzzii: 0.023 from the rounded table, 0.006 in the text).
