# Malaria vector software — tool reference

Short reference for the three cloud-backed / data-analysis Python packages this project depends on for target selection and downstream analysis. Detailed agent-facing usage notes for each tool live under `skills/<tool>/` in this repository; this file is the one-page overview.

---

## malariagen_data

**Programmatic access to MalariaGEN genomic surveillance data.**

- Python package: `pip install malariagen_data`.
- Streams data directly from Google Cloud Storage (Zarr-backed); no local files required.
- Covers Anopheles vector releases — **Ag3** (*An. gambiae s.l.*), **Af1** (*An. funestus*), **Amin1**, **Adir1** — and Plasmodium parasite releases (Pf7/Pf8/Pf9, Pv4).
- Entry point for this project is the `Ag3` class: SNP genotypes + allele frequencies, CNVs, haplotypes, sample metadata, AIMs, karyotypes, reference genome + gene features, selection scans, diversity statistics.
- Authentication: `gcloud auth application-default login` outside Colab.

Used by this project for: SNP allele-frequency queries (primer binding-site filtering), AIM selection (Module E), variable-locus selection (Module F), and reference/GFF ingest.

Skill reference: `skills/malariagen_data/` — see `index.md`, `snp_data.md`, `frequencies.md`, `aim_data.md`, `selection_scans.md`.

---

## selection-atlas

**Catalogue of positive selection signals across the major malaria vectors.**

- Authors: Nagi, Harding, Lawniczak, Donnelly, Miles.
- Covers the *An. gambiae* complex and *An. funestus*, built from whole-genome Ag3 / Af1 data.
- Runs H12, G123 and iHS selection scans per cohort; a peak-finding algorithm applied to H12 produces a single flat signals table.
- **No dedicated Python API.** The primary entry point is the signals CSV: `h12-signal-detection-all.csv`. Load with pandas. Deeper follow-up uses `malariagen_data` directly.
- Each row = one H12 signal in one cohort on one contig.

Used by this project for: Module B target selection — top sweep windows flagging resistance-associated loci.

Skill reference: `skills/selection-atlas/` — see `index.md`, `signals.md` (column schema), `queries.md`.

---

## AnoExpress

**Meta-analysis of Anopheles gene expression in insecticide-resistance studies.**

- Python package: `pip install anoexpress`, then `import anoexpress as ae`.
- Aggregates RNA-Seq (and microarray via IR-Tex) studies of insecticide resistance in *An. gambiae s.l.* and *An. funestus*.
- Expression data streamed remotely from GitHub; no local files required.
- Key API parameter: `analysis` — selects species set (`gamb_colu`, `gamb_colu_arab`, `fun`, etc.). Analyses with more species contain fewer genes (ortholog constraint).
- Functions cover: differential expression tables, candidate-gene lookup, enrichment, plotting.

Used by this project for: Module C target selection — consistently DE genes (resistant vs susceptible) flagging new resistance candidates.

Skill reference: `skills/AnoExpress/` — see `index.md`, `data.md`, `candidates.md`, `enrichment.md`.

---

## How they fit together in this project

| Module | Purpose                           | Source                                                      |
|--------|-----------------------------------|-------------------------------------------------------------|
| A      | Known resistance loci             | Literature + AnoPrimer                                      |
| B      | Sweep-derived resistance          | **selection-atlas** (CSV → pandas)                          |
| C      | Expression-derived resistance     | **AnoExpress** (Python API)                                 |
| D      | Orthologue-transferred resistance | irtho + discovery-species inputs                            |
| E      | AIMs                              | **malariagen_data** (`Ag3`, allele-freq contrasts)          |
| F      | Variable loci                     | **malariagen_data** (`Ag3`, π / diversity)                  |

All three tools ultimately draw on MalariaGEN Ag3 / Af1 cloud releases; selection-atlas and AnoExpress are curated views over that underlying data.
