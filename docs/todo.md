# Bouaké paper — remaining analyses

Phased plan for final analyses prior to submission. Each phase has a dedicated plan file.

## Phase A — Vgsc F1529C (demoted to a short section)
Decision: keep F1529C in the Bouaké paper as one paragraph in "Variants of interest" plus
Supplementary Figure 2 (Bouaké reads and haplotype background). The Ag3 survey (maps, time course,
country breakdown, haplotype sharing statistics) is saved for a separate paper; its code and results are
kept (`analysis/f1529c/`, `results/f1529c/`, `figures_ms/figure_f1529c.py`, `figure_f1529c_cohorts.py`).

- [x] Paragraph, Methods paragraph and Supplementary Figure 2 written into the manuscript
- [x] All Ag3 calls restricted to sample sets with `unrestricted_use` true (`src/bouake/ag3_sets.py`)
- [ ] Replace the `(ref)` placeholders (pyrethroid resistance of Aedes F1534C; scaffold-phasing method)
- [ ] Framing: Abstract, Discussion and title still need a sentence on why a pyrethroid-target mutation
      appears in a pirimiphos-methyl paper (multiple resistance in the study population)
- [ ] Separate F1529C paper: survey, time course, single-origin test, a pyrethroid association if available

## Phase B — ABCs (gambiae) + OBPs (coluzzii) deep dive
Plan: `docs/phase-b-abc-obp-plan.md`

Refocused mid-phase onto a single candidate gene, *ABCH1*, after the family-level
analyses came back null (see below). Figure 5 and its Results section are written.

- [x] B1. Extract DE ABCs (gambiae) + DE OBPs (coluzzii) from existing DE tables
- [x] B2. AnoExpress cross-study context for ABC family (gambiae focus)
- [x] B3. AnoExpress cross-study context for OBP family (coluzzii focus)
- [x] B4. Cross-study reproducibility — done as T4 (Mann–Whitney rank shift), all null
- [x] B5. Selection-atlas sweep check at top ABC and OBP loci
- [x] B6. Figure 5 — *ABCH1* expression (A), H12 sweep (B), locus zoom (C)
- [x] B8. Real H12 scan for 8 Ag3 coluzzii cohorts (run on hyperion; `results/abc_obp/h12/`)
- [ ] B7. Supplementary tables — ranked DE + cross-study support
- [x] B9. (done) **Fix the ABC family definition** — `GO:0042626` pulls in non-ABC ATPases
      (AGAP002858 Na+/K+-ATPase, AGAP000523 F-type H+-ATPase and 8 others are counted
      among the "30 sig-up ABCs"). Restrict to PFAM PF00005 + `Abc*` names and
      regenerate `results/abc_obp/bouake_de_abc_*.csv`. Does not affect Figure 5.
- [x] B10. (done) Decide whether the Figure 2 GSEA paragraph needs revising — OBPs are up in
      both species (17 gambiae / 16 coluzzii, sign concordance 0.94), so the
      "gambiae = ABC, coluzzii = OBP" dichotomy does not hold as a between-species
      claim. The ABC direction reversal (30 up/5 dn vs 6 up/26 dn) does hold.

### Phase B null results (reported, not hidden)
- Set-level sweep enrichment among sig-up ABCs: **P = 0.18** after collapsing tandem
  arrays (the uncollapsed gene-level P = 0.019 was inflated by one 4-gene ABCC array).
- Genomic clustering of DE response (T3-redux): null for every family/species.
  The earlier coluzzii-ABC P = 0.0001 was a 2-gene-adjacency artefact; now P = 1.00.
- Cross-study rank shift (T4): null for all six family × species sets, padj > 0.98.
- OBP selection: 0 of 75 OBPs overlap any H12 peak focus in West Africa.

## Infrastructure
- [x] Pixi init at repo root with required deps (zarr pinned <3 for malariagen_data compatibility)
- [x] Confirm `gcloud auth application-default login` available (needed by malariagen_data)
- [x] Notebook + results dirs under `analysis/f1529c/`, `analysis/abc_obp/`, `results/f1529c/`, `results/abc_obp/`
- Env entry point: `pixi run python ...` or `pixi run lab` for Jupyter

## ABC family fix (B9/B10) — done 2026-10-07
Strict family = PF00005 or `Abc*` name (56 genes; was 99 with GO:0042626 and four TM-domain PFAMs).
Setup and every dependent script now use it (`analysis/abc_obp/setup.py`; caches live in
`results/abc_obp/cache/`; `b9_abc_tables.py` regenerates the ABC tables). Corrected results:
- gambiae 24 sig ABCs: 20 up / 4 down; coluzzii 12 sig: 5 up / 7 down (Fisher P = 0.02)
- the 20 up-regulated gambiae ABCs are clustered (median NN 19 kb vs 1.4 Mb; adj P = 0.008), two tandem arrays on 3R
- set-level sweep overlap, tandem arrays collapsed: 12 loci, 3 swept, P = 0.27, OR 1.7 (null; was P = 0.18, OR 2.07)
- ABCH1 is the only ABC gene with a focus-level sweep signal; cross-study rank shift (T4) still null
- the OBP direction result is unchanged: OBPs up in both species (17 / 16)
Removed (they used the old ABC set): `figures_ms/b1_*abc*`, `b2_abc*`, `b5_abc*` (png and html) and
`results/abc_obp/bouake_abc_gambiae_vs_anoexpress.csv`. Notebooks 01-03, 05 and 06 in `analysis/abc_obp/` still
contain the old definition and would regenerate those files; the maintained path is `setup.py` and the scripts.
Also corrected in the paper: the phthalate-binding "padj" quoted the raw P (padj is 1.1e-3).

## Embargo fix (Ag3 sample sets) — done 2026-10-07
All Ag3 calls in `analysis/abc_obp/` now pass `sample_sets=unrestricted_sample_sets(...)` (`src/bouake/ag3_sets.py`,
moved from `analysis/f1529c/`). The 55 selection-atlas cohorts contain no restricted samples, so the H12 signals,
sweep overlaps and the 8-versus-0 cohort count were already clean. Rerun on unrestricted sets only: haplotype clustering
(5,804 haplotypes; swept clade 663, all coluzzii), PBS, H1X, copy number (3,991 of 3,997 at two copies), Figure 4.
Outputs now write straight to `results/abc_obp/<name>/`. Removed the unused `njt_haplotypes.py`/`njt_only.py` and the
broken NJT block of the CNV script (`cnv_njt_h1x.py` renamed `cnv_h1x.py`). `h12_abch1.py` and `extra_stats.py` were not
rerun (their cohorts are unrestricted) and still write to `./h12_out`, `./stats_out`.
