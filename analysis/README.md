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
