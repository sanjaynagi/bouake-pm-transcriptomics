"""Additional selection statistics across 2R for the Figure 5 multi-track panel.

H12 is already cached; this adds iHS and G123 for the same An. coluzzii cohorts, so the
locus can be shown the way Ingham et al. (2020) show SAP2 — several statistics stacked over
a shared genomic axis.

Run on hyperion. Writes to ./stats_out/.
"""

import os
import numpy as np
import malariagen_data

CONTIG = "2R"
ANALYSIS = "gamb_colu"
OUT = "stats_out"

COHORTS = [
    ("BF-02_Comoe_colu_2012",              "year"),
    ("BF-09_Houet_colu_2012_Q3",           "quarter"),
    ("BF-02_Comoe_colu_2015",              "year"),
    ("BF-02_Comoe_colu_2016",              "year"),
    ("BF-09_Houet_colu_2014_Q3",           "quarter"),
    ("BF-02_Comoe_colu_2011",              "year"),
    ("ML-3_Yanfolila_colu_2012_Q4",        "quarter"),
    ("GH-AH_Adansi-Akrofuom_colu_2018_Q4", "quarter"),
]

os.makedirs(OUT, exist_ok=True)
ag3 = malariagen_data.Ag3(results_cache="~/ag3_cache_h12", pre=True)

for cohort_id, level in COHORTS:
    col = "cohort_admin2_quarter" if level == "quarter" else "cohort_admin2_year"
    query = f"{col} == '{cohort_id}'"

    dest = f"{OUT}/{cohort_id}.ihs.npz"
    if not os.path.exists(dest):
        try:
            print(f"[ihs ] {cohort_id}", flush=True)
            x, ihs = ag3.ihs_gwss(contig=CONTIG, analysis=ANALYSIS, sample_query=query,
                                  window_size=200, min_cohort_size=15, max_cohort_size=50,
                                  random_seed=42)
            np.savez_compressed(dest, x=x, ihs=ihs, cohort=cohort_id)
            print(f"[ok  ] ihs {cohort_id}: {len(x)} windows", flush=True)
        except Exception as e:
            print(f"[FAIL] ihs {cohort_id}: {type(e).__name__}: {e}", flush=True)

    dest = f"{OUT}/{cohort_id}.g123.npz"
    if not os.path.exists(dest):
        try:
            print(f"[g123] {cohort_id}", flush=True)
            x, g123 = ag3.g123_gwss(contig=CONTIG, window_size=1000, sites=ANALYSIS,
                                    site_mask="gamb_colu", sample_query=query,
                                    min_cohort_size=15, max_cohort_size=50, random_seed=42)
            np.savez_compressed(dest, x=x, g123=g123, cohort=cohort_id)
            print(f"[ok  ] g123 {cohort_id}: {len(x)} windows", flush=True)
        except Exception as e:
            print(f"[FAIL] g123 {cohort_id}: {type(e).__name__}: {e}", flush=True)

print("done", flush=True)
