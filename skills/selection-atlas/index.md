---
name: selection-atlas
description: Overview and index for the Malaria Vector Selection Atlas — a catalogue of positive selection signals across the Anopheles gambiae complex and An. funestus. Use this skill to understand the resource, the signal-detection approach, and to navigate to the right reference page for interrogating the signal table and running follow-up analyses with malariagen_data.
---

# Selection Atlas

**An atlas of positive selection in the genomes of major malaria vectors**
(Nagi, Harding, Lawniczak, Donnelly, Miles)

A catalogue of recent positive selection in natural populations of the major African malaria vectors (*An. gambiae* complex; *An. funestus*), built from whole-genome sequence data held in the MalariaGEN *Ag3* and *Af1* data releases. Whole-genome H12, G123 and iHS selection scans are run for every cohort, and a peak-finding algorithm is applied to H12 to produce a *single flat table* of selection signals. That table is the primary entry point for programmatic use.

There is **no dedicated Python API** — all interrogation is done by loading the signals CSV into pandas and, where deeper analysis is wanted, using `malariagen_data` (`Ag3` / `Af1`) to re-run scans, pull haplotype data, fetch genome features, etc.

---

## The signals table

**File:** `h12-signal-detection-all.csv` (at the repo root).

Each row is one H12 selection signal fitted in one cohort on one contig. The table is the deduplicated union of per-cohort signal detection outputs (see `workflow/gwss/notebooks/h12-signal-detection.ipynb`).

See [signals.md](signals.md) for the full column schema.

---

## Key concepts

### How signals are detected (in brief)

For each cohort × contig, H12 is computed genome-wide (`malariagen_data.Ag3.h12_gwss`) with a cohort-calibrated window size, then:

1. **Hampel filter** removes outlier windows.
2. Physical positions (bp) are mapped to **genetic positions (cM)** via a recombination map.
3. A **skewed exponential-decay peak** is fit in 1 cM sliding windows alongside a null constant model.
4. A signal is emitted where **ΔAIC = AIC(null) − AIC(peak) > 1000** *and* peak H12 > 0.1.
5. Overlapping signals are deduplicated, keeping the one with the largest ΔAIC.

The fitted peak gives both a location (`gcenter`/`pcenter`) and a *width* via the decay constant, from which three concentric intervals are derived (see below).

### Signal intervals (focus / span1 / span2)

Each signal has three concentric intervals around the fitted peak center, reported in both genetic (`g*`, cM) and physical (`p*`, bp) coordinates:

| Interval | Width | Use for |
|---|---|---|
| `focus`  | center ± 0.25 × decay | tightest core — most likely to contain the selected locus |
| `span1`  | center ± 1 × decay    | conservative signal extent |
| `span2`  | center ± 2 × decay    | broad footprint — useful for gene overlap / LD |

The peak is **skewed**, so `decay_left` and `decay_right` can differ. All three intervals respect that skew.

### `cohort_id` schema

Cohort identifiers encode country, admin unit, taxon, year and (usually) quarter:

```
{country}-{admin1_iso}_{admin2}_{taxon}_{year}[_Q{quarter}]
```

- `country`    — ISO alpha-2 code (e.g. `BF`, `CI`, `ML`, `UG`, `TZ`)
- `admin1_iso` — admin-1 ISO code within the country (e.g. `BF-09`)
- `admin2`    — admin-2 name, hyphenated (e.g. `Houet`, `Adansi-South`)
- `taxon`    — `gamb` (*gambiae*), `colu` (*coluzzii*), `arab` (*arabiensis*), `biss` (*bissau*), or funestus equivalents
- `year`     — 4-digit year of collection
- `quarter`  — optional `Q1`–`Q4`; absent if samples were not provided with a quarter

Examples: `BF-09_Houet_colu_2012_Q3`, `ML-2_Kangaba_gamb_2004_Q3`, `GW-BS_Bissau-Autonomous-Sector_biss_2010`.

The `cohort_id` is also a valid value to query against `cohort_admin2_quarter` / `cohort_admin2_year` columns in `malariagen_data` sample metadata — see [followup.md](followup.md).

### Contigs

Cohorts in the *Ag3* atlas are scanned on the compound arms `2RL`, `3RL` and `X`. Cohorts in the *Af1* atlas use `2RL`, `3RL` and `X` with the *An. funestus* reference.

---

## Skills index

| Skill file | Contents |
|---|---|
| [signals.md](signals.md)  | Full column schema for `h12-signal-detection-all.csv` and how to load it |
| [queries.md](queries.md)  | Recipes for common questions — genes under sweeps, strongest signals, parallel selection across cohorts |
| [followup.md](followup.md) | Running deeper analyses on a signal with `malariagen_data` (re-run H12, haplotype networks, gene features, SNP frequencies) |

---

## Companion skills

For follow-up analysis, these skills (installed separately) are the main toolkit:

- **malariagen_data** (`skills-malariagen_data/`) — genome-wide selection scans, haplotype data, SNP calls, sample metadata, frequencies, gene features.
- **AnoExpress** (`AnoExpress/skills/`) — cross-referencing selection signals with insecticide-resistance differential-expression data.

---

## Citation

Nagi SC, Harding NJ, Lawniczak MKN, Donnelly MJ, Miles A. *An atlas of positive selection in the genomes of major malaria vectors.* See `selection-atlas ms.md` in this repo for the manuscript.
