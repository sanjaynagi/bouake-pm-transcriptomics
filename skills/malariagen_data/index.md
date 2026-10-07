---
name: malariagen-data
description: Overview and index for the malariagen_data Python package — genomic data access for Anopheles mosquitoes (Ag3, Af1, Amin1, Adir1) and Plasmodium parasites (Pf7/Pf8/Pf9, Pv4). Use this skill to understand what the package does and navigate to the right API reference.
---

# malariagen_data

**Genomic surveillance data for malaria vectors and parasites**

A Python package providing programmatic access to data from the Malaria Genomic Epidemiology Network (MalariaGEN). All data is streamed directly from Google Cloud Storage — no local data files are required.

## Installation

```bash
pip install malariagen_data
```

### Authentication

From Google Colab, authentication is initiated automatically.

Outside Colab, install the Google Cloud CLI and run:

```bash
gcloud auth application-default login
```

## Entry points

| Class | Species | Example |
|---|---|---|
| `Ag3` | *Anopheles gambiae* complex | `ag3 = malariagen_data.Ag3()` |
| `Af1` | *Anopheles funestus* subgroup | `af1 = malariagen_data.Af1()` |
| `Amin1` | *Anopheles minimus* | `amin1 = malariagen_data.Amin1()` |
| `Adir1` | *Anopheles dirus* complex | `adir1 = malariagen_data.Adir1()` |
| `Pf7` / `Pf8` / `Pf9` | *Plasmodium falciparum* | `pf7 = malariagen_data.Pf7()` |
| `Pv4` | *Plasmodium vivax* | `pv4 = malariagen_data.Pv4()` |

## Common parameters

### Sample selection

| Parameter | Type | Description |
|---|---|---|
| `sample_sets` | `str` or `list[str]` | One or more sample set IDs (e.g. `"AG1000G-AO"`) or a release tag (e.g. `"3.0"`) |
| `sample_query` | `str` | Pandas query string applied to sample metadata (e.g. `"country == 'Ghana'"`) |
| `sample_indices` | `list[int]` | Integer positional indices into the sample array — mutually exclusive with `sample_query` |
| `cohort_size` | `int` | Exact number of samples drawn at random per cohort |
| `min_cohort_size` / `max_cohort_size` | `int` | Bounds on random subsampling |

### Genomic regions

Regions can be specified as:
- Contig string: `"2L"`
- Region string: `"2L:1000000-2000000"`
- `Region` named tuple

### Site masks

Site masks filter to high-quality accessible sites. Common values for Ag3: `"gamb_colu"`, `"gamb_colu_arab"`, `"arab"`. Default is `None` (no filtering) for most functions; many functions default to `base_params.DEFAULT` which selects the species-appropriate mask.

### Data format parameters

| Parameter | Type | Description |
|---|---|---|
| `inline_array` | `bool` | If `True`, load arrays inline (faster for small data) |
| `chunks` | `str` or `tuple` | Dask chunking strategy — `"native"` uses stored chunk sizes |

---

## Skills index

| Skill file | Contents |
|---|---|
| [setup.md](setup.md) | Initialising API objects and basic data access |
| [sample_metadata.md](sample_metadata.md) | Sample metadata, cohort tables, maps |
| [reference_genome.md](reference_genome.md) | Genome sequence, gene features, transcript plots |
| [snp_data.md](snp_data.md) | SNP calls, allele counts, site annotations |
| [haplotypes.md](haplotypes.md) | Phased haplotype data |
| [aim_data.md](aim_data.md) | Ancestry-informative markers — AIM calls and heatmap (Ag3) |
| [karyotype.md](karyotype.md) | Chromosomal inversion karyotyping (Ag3) |
| [cnv_data.md](cnv_data.md) | CNV HMM, coverage calls, gene copy number |
| [frequencies.md](frequencies.md) | SNP, amino-acid, CNV and haplotype frequencies |
| [selection_scans.md](selection_scans.md) | H12, H1X, G123, IHS, XPEHH genome-wide scans |
| [population_structure.md](population_structure.md) | PCA, NJT, FST, diplotype and haplotype clustering |
| [diversity.md](diversity.md) | Diversity stats, heterozygosity, runs of homozygosity |
| [phenotypes.md](phenotypes.md) | Phenotype data, IGV browser, cross metadata (Ag3) |
