"""Compute H12 across 2R for the An. coluzzii cohorts carrying a selection-atlas signal that
spans ABCH1 (AGAP002638, 2R:24,825,285-24,848,323), for Figure 5B.

Run on the LSTM cluster (hyperion), where gcloud application-default credentials are available.
Writes one npz per cohort into ./h12_out/.
"""

import os
import numpy as np
import malariagen_data

WINDOW_SIZE = 1000          # consistent across cohorts so traces are directly comparable
CONTIG = "2R"
ANALYSIS = "gamb_colu"      # phasing analysis covering gambiae/coluzzii

# Cohorts whose span1 interval overlaps ABCH1 (from h12-signal-detection-all.csv).
# focus=True marks the five whose peak FOCUS also sits on the gene.
COHORTS = [
    ("BF-02_Comoe_colu_2012",              "year",    True),
    ("BF-09_Houet_colu_2012_Q3",           "quarter", True),
    ("BF-02_Comoe_colu_2015",              "year",    True),
    ("BF-02_Comoe_colu_2016",              "year",    True),
    ("BF-09_Houet_colu_2014_Q3",           "quarter", True),
    ("BF-02_Comoe_colu_2011",              "year",    False),
    ("ML-3_Yanfolila_colu_2012_Q4",        "quarter", False),
    ("GH-AH_Adansi-Akrofuom_colu_2018_Q4", "quarter", False),
]

OUT = "h12_out"
os.makedirs(OUT, exist_ok=True)

ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)

for cohort_id, level, is_focus in COHORTS:
    dest = f"{OUT}/{cohort_id}.npz"
    if os.path.exists(dest):
        print(f"[skip] {cohort_id}", flush=True)
        continue
    col = "cohort_admin2_quarter" if level == "quarter" else "cohort_admin2_year"
    query = f"{col} == '{cohort_id}'"
    print(f"[run ] {cohort_id} ({query})", flush=True)
    try:
        x, h12, contigs = ag3.h12_gwss(
            contig=CONTIG,
            window_size=WINDOW_SIZE,
            analysis=ANALYSIS,
            sample_query=query,
            min_cohort_size=15,
            max_cohort_size=50,
            random_seed=42,
        )
        np.savez_compressed(dest, x=x, h12=h12, cohort=cohort_id, focus=is_focus,
                            window_size=WINDOW_SIZE)
        print(f"[ok  ] {cohort_id}: {len(x)} windows, max H12 {np.nanmax(h12):.3f}", flush=True)
    except Exception as e:
        print(f"[FAIL] {cohort_id}: {type(e).__name__}: {e}", flush=True)

print("done", flush=True)
