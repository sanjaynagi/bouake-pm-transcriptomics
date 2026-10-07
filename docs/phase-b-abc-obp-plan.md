# Phase B — ABCs (gambiae) + OBPs (coluzzii) plan

## Problem

The GSEA results (Figure 2) show a strikingly species-specific pattern in Bouaké:

- **An. gambiae** — ABC transporters / transmembrane efflux dominate (`ABC_tran` PFAM padj = 1.3×10⁻⁹; GO:0055085 padj = 3.5×10⁻¹²; KEGG aga02010 ABC transporters padj = 1.5×10⁻³).
- **An. coluzzii** — chemosensory / odorant-binding proteins dominate (`PBP_GOBP` PFAM padj = 6.9×10⁻⁸; GO:0005549 odorant binding padj = 7.0×10⁻⁷; `DUF753` padj = 9.5×10⁻¹⁶).

The manuscript currently reports these enrichments in one paragraph each, but Figure 5 ("ORs and ABCs") is a blank placeholder. Three open questions:

1. **Which specific ABCs and OBPs drive the enrichment in Bouaké?** Not yet listed individually.
2. **Are these signals reproducible in the wider Anopheles resistance literature, or novel to Bouaké?** AnoExpress is the direct tool for this.
3. **Do the ABC / OBP loci show evidence of positive selection in wild An. gambiae s.l. populations?** Selection-atlas is the direct tool.

## Proposed solution

### B1. Bouaké DE extraction

Load `results/BouakePM_diffexp.xlsx` (and `BouakePM_diffexp_P450_ABC_GST_COE.xlsx`). For each species-specific contrast (coluzziiCont vs coluzziiPM; gambiaeCont vs gambiaePM):

- Identify DE ABCs using PFAM `ABC_tran` + GO `GO:0042626` / `GO:0055085` via `ae.load_annotations()`.
- Identify DE OBPs using PFAM `PBP_GOBP` + GO `GO:0005549`.
- Export ranked tables (log2FC, padj, annotation) into `results/abc_obp/`.

### B2. AnoExpress cross-study context — ABCs (gambiae)

```python
ae.plot_gene_family_expression(
    gene_identifier="ABC_tran",
    analysis="gamb_colu_arab_fun",
    title="ABC transporters across Anopheles IR studies",
    plot_type="strip",
    sort_by="median",
)
```

Then `ae.load_candidates(analysis="gamb_colu_arab_fun", query_annotation="ABC_tran", func=np.nanmedian)` to rank ABCs by median fold change across the meta-analysis. Cross-reference Bouaké gambiae hits against this ranking → flag each Bouaké ABC as `shared_candidate` or `novel`.

### B3. AnoExpress cross-study context — OBPs (coluzzii)

Same approach with `gene_identifier="PBP_GOBP"`. Note the OBP story is more surprising — OBPs are not canonical resistance genes, so a strong Bouaké signal that is absent in the AnoExpress meta-analysis would itself be a headline.

Also run with `gene_identifier="GO:0005549"` (odorant binding) as a sensitivity check.

### B4. Consistent-genes overlap

```python
ae.consistent_genes(analysis="gamb_colu_arab_fun", direction="up", n_experiments=5)
```

Intersect the resulting gene set with Bouaké DE ABCs and DE OBPs. Genes present in both = robustly reproducible resistance candidates. Genes Bouaké-only = novel / site-specific signal.

### B5. Selection-atlas sweep check

Load `skills/selection-atlas/h12-signal-detection-all.csv`. For each top DE ABC / OBP, lookup the gene coordinates via `ag3.geneset()` then query:

```python
signals.query(
    "contig == @contig and span1_pstart <= @gstop and span1_pstop >= @gstart"
).sort_values("delta_i", ascending=False)
```

Filter cohorts to West Africa (`country in ['CI','BF','GH','ML','TG','BJ','BN','NG']`). Expected values: most ABCs in efflux / detox show modest sweeps; OBP 2L cluster may show no sweep (would align with expression-only response rather than hard selection).

Tag each DE gene with: `under_sweep_in_west_africa` (Y/N), and if yes — cohorts, delta_i, stat_max.

### B6. Figure 5 composition

Two-panel figure:

- **5A (gambiae ABC panel)** — heatmap of top DE ABCs in Bouaké gambiae (rows) with AnoExpress meta-analysis fold-change distribution (boxplot / strip next to each row).
- **5B (coluzzii OBP panel)** — same for top DE OBPs in coluzzii.
- Both annotated with a side column indicating sweep overlap (black/grey ticks).

### B7. Supplementary tables

- `Supp_Table_ABC_DE.xlsx` — Bouaké gambiae DE ABCs with log2FC, padj, AnoExpress median FC, n_AnoExpress_sig, sweep overlap.
- `Supp_Table_OBP_DE.xlsx` — same structure for coluzzii OBPs.

## Risks and tradeoffs

- **PFAM ABC_tran is broad** — includes non-efflux ABCs (e.g. transporters unrelated to xenobiotics). May want to restrict to ABCB / ABCG / ABCH subfamilies, which are the xenobiotic-relevant ones (Ingham et al. 2021). Use GeneName patterns (e.g. `Abcb`, `Abcg`, `Abch`).
- **OBP genomic clustering** — coluzzii DE OBPs are likely in the 2L chemosensory array; if so, enrichment could reflect one trans-acting regulator rather than many independent signals. Worth noting in the Discussion, not a blocker.
- **AnoExpress coverage** — meta-analysis is dominated by pyrethroid resistance studies; absence of a Bouaké signal in AnoExpress may reflect insecticide specificity (PM-only) rather than true novelty. Flag this in the output.
- **Selection-atlas coverage** — Ag3 data; Bouaké 2019 samples are not in Ag3. Sweep analysis is therefore "are these loci under selection in West African wild populations generally?", not "in Bouaké specifically".

## Success criteria

- [ ] Ranked DE tables for ABCs (gambiae) and OBPs (coluzzii) with annotation + padj + log2FC.
- [ ] AnoExpress cross-study fold-change plot for each family.
- [ ] Consistent-genes overlap table — how many Bouaké hits are shared vs novel.
- [ ] Selection-atlas sweep overlap table for each DE gene.
- [ ] Figure 5 completed (two-panel) replacing the blank placeholder.
- [ ] Two supplementary tables.
- [ ] Draft paragraphs for Results extending the current GSEA section, and for Discussion.

## Deliverables

- `notebooks/abc_obp/01_bouake_de_extraction.ipynb`
- `notebooks/abc_obp/02_anoexpress_abc.ipynb`
- `notebooks/abc_obp/03_anoexpress_obp.ipynb`
- `notebooks/abc_obp/04_sweep_overlap.ipynb`
- `notebooks/abc_obp/05_figure5.ipynb`
- `results/abc_obp/bouake_de_abc.csv`
- `results/abc_obp/bouake_de_obp.csv`
- `results/abc_obp/sweep_overlap.csv`
- `figures_ms/figure5_abc_obp.{svg,png}`
- `bouake_supplement/Supp_Table_ABC_DE.xlsx`
- `bouake_supplement/Supp_Table_OBP_DE.xlsx`
