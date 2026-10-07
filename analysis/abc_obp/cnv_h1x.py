"""Two additions to the ABCH1 locus analysis.

1. CNV — gene amplification is the recurrent mechanism at Anopheles detox/transport loci, so
   a CNV under the sweep summit would be a strong candidate causal variant. Copy-number state
   per sample, plus CNV frequency by admin2 x year to see whether it rises through time.
2. H1X — cross-population haplotype homozygosity, An. coluzzii vs An. gambiae. ABCH1 is
   upregulated in both species in Bouaké, so this tests whether the swept haplotype is shared
   between them (adaptive introgression) or has arisen independently.

Writes to results/abc_obp/cnv/.
"""

import os
import numpy as np
import pandas as pd
import malariagen_data
from pathlib import Path

from bouake.ag3_sets import unrestricted_sample_sets  # noqa: E402

RESULTS = Path(__file__).resolve().parents[2] / "results" / "abc_obp"

CONTIG = "2R"
ANALYSIS = "gamb_colu"
LOCUS = "2R:24,780,000-24,900,000"     # ABCH1 through Or38 plus flanks
ABCH1_REGION = "2R:24,825,285-24,848,323"
OUT = str(RESULTS / "cnv")
WA = ["Burkina Faso", "Mali", "Ghana", "Cote d'Ivoire"]

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)
SAMPLE_SETS = unrestricted_sample_sets(ag3=ag3, taxa=("coluzzii", "gambiae"))
colu_wa = "taxon == 'coluzzii' and country in {!r}".format(WA)

# ── 1. CNV ──────────────────────────────────────────────────────────────
try:
    print("[cnv] gene_cnv over the locus", flush=True)
    ds = ag3.gene_cnv(region=LOCUS, sample_query=colu_wa, sample_sets=SAMPLE_SETS, max_coverage_variance=0.2)
    print("  dims:", dict(ds.sizes), flush=True)
    genes = ds["gene_id"].values
    cn = ds["CN_mode"].values                      # genes x samples
    meta = pd.DataFrame(dict(sample_id=ds["sample_id"].values,
                             country=ds["sample_country"].values if "sample_country" in ds else None))
    amp = pd.DataFrame(cn, index=genes, columns=ds["sample_id"].values)
    amp.to_csv(f"{OUT}/gene_cn_mode.csv")
    frac_amp = (amp > 2).mean(axis=1).sort_values(ascending=False)
    print("\n  fraction of samples with CN > 2, top genes:", flush=True)
    print(frac_amp.head(12).to_string(), flush=True)
    print(f"\n  ABCH1 (AGAP002638): "
          f"{frac_amp.get('AGAP002638', float('nan')):.4f} amplified", flush=True)
except Exception as e:
    print(f"[FAIL] gene_cnv: {type(e).__name__}: {e}", flush=True)

# CNV frequency through time — admin2 x year cohorts
for cohort_key in ("admin2_year", "admin1_year"):
    try:
        print(f"\n[cnv] gene_cnv_frequencies by {cohort_key}", flush=True)
        fr = ag3.gene_cnv_frequencies(region=LOCUS, cohorts=cohort_key,
                                      sample_query=colu_wa, sample_sets=SAMPLE_SETS, min_cohort_size=10,
                                      drop_invariant=False, include_counts=True)
        fr.to_csv(f"{OUT}/gene_cnv_freq_{cohort_key}.csv")
        abch1 = fr[fr.index.get_level_values(0).astype(str).str.contains("AGAP002638")] \
            if fr.index.nlevels > 1 else fr[fr.index.astype(str).str.contains("AGAP002638")]
        print(abch1.to_string()[:3000], flush=True)
        break
    except Exception as e:
        print(f"[FAIL] gene_cnv_frequencies {cohort_key}: {type(e).__name__}: {e}", flush=True)

# ── 2. H1X, coluzzii vs gambiae in Burkina Faso ─────────────────────────
try:
    print("\n[h1x] coluzzii vs gambiae, Burkina Faso", flush=True)
    x, h1x, contigs = ag3.h1x_gwss(
        contig=CONTIG, window_size=1000, analysis=ANALYSIS,
        cohort1_query="taxon == 'coluzzii' and country == 'Burkina Faso'",
        cohort2_query="taxon == 'gambiae' and country == 'Burkina Faso'",
        sample_sets=SAMPLE_SETS,
        min_cohort_size=15, max_cohort_size=50, random_seed=42)
    np.savez_compressed(f"{OUT}/h1x_BF_colu_vs_gamb.npz", x=x, h1x=h1x)
    m = (x > 24.78e6) & (x < 24.90e6)
    print(f"  windows {len(x)}; max H1X at locus {np.nanmax(h1x[m]):.3f}; "
          f"genome median {np.nanmedian(h1x):.3f}", flush=True)
except Exception as e:
    print(f"[FAIL] h1x: {type(e).__name__}: {e}", flush=True)

print("\ndone", flush=True)
