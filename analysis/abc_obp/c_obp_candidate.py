"""C. Is there an OBP equivalent of ABCH1 — strongly DE in Bouaké AND under a replicated,
focus-level sweep in wild West African populations?

Two levels of evidence, matching how ABCH1 qualifies:
  span-level  — gene falls inside the broad signal interval (weak; ABCC array only qualifies here)
  focus-level — gene overlaps the signal FOCUS, i.e. the inferred sweep centre (strong; ABCH1)
"""

import numpy as np
import pandas as pd
from setup import gam, col, FAM, genes_df, sig_wa

pd.set_option("display.width", 260)
pd.set_option("display.max_colwidth", 60)

coords = genes_df.set_index("GeneID")


def signal_hits(gene_id, level="span"):
    """Signals overlapping a gene, at span or focus resolution."""
    if gene_id not in coords.index:
        return sig_wa.iloc[0:0]
    g = coords.loc[gene_id]
    s = sig_wa[sig_wa["contig"] == g["atlas_contig"]]
    lo, hi = (("span1_pstart", "span1_pstop") if level == "span"
              else ("focus_pstart", "focus_pstop"))
    return s[(s[lo] <= g["atlas_end"]) & (s[hi] >= g["atlas_start"])]


def profile(gene_ids, family_label):
    rows = []
    for gid in gene_ids:
        g = gam[gam["GeneID"] == gid]
        c = col[col["GeneID"] == gid]
        span, focus = signal_hits(gid, "span"), signal_hits(gid, "focus")
        name = (g["GeneName"].iloc[0] if len(g) and pd.notna(g["GeneName"].iloc[0])
                else (coords.loc[gid, "Name"] if gid in coords.index else np.nan))
        rows.append(dict(
            family=family_label, GeneID=gid, name=name if pd.notna(name) else gid,
            contig=coords.loc[gid, "contig"] if gid in coords.index else None,
            l2fc_gam=g["log2FoldChange"].iloc[0] if len(g) else np.nan,
            padj_gam=g["padj"].iloc[0] if len(g) else np.nan,
            l2fc_col=c["log2FoldChange"].iloc[0] if len(c) else np.nan,
            padj_col=c["padj"].iloc[0] if len(c) else np.nan,
            n_span=len(span), n_focus=len(focus),
            focus_cohorts=len(focus["cohort_id"].unique()),
            focus_countries=",".join(sorted(focus["country"].unique())),
            max_delta_i=focus["delta_i"].max() if len(focus) else np.nan,
        ))
    return pd.DataFrame(rows)


obp = profile(sorted(FAM["OBP"]), "OBP")
abc = profile(sorted(FAM["ABC"]), "ABC")

obp["sig_both"] = (obp["padj_gam"] < 0.05) & (obp["padj_col"] < 0.05) & \
                  (obp["l2fc_gam"] > 0) & (obp["l2fc_col"] > 0)
abc["sig_both"] = (abc["padj_gam"] < 0.05) & (abc["padj_col"] < 0.05) & \
                  (abc["l2fc_gam"] > 0) & (abc["l2fc_col"] > 0)

print("=" * 120)
print("C1. OBP family — every gene with ANY sweep overlap (span or focus)")
print("=" * 120)
hit = obp[(obp["n_span"] > 0) | (obp["n_focus"] > 0)].sort_values("n_focus", ascending=False)
print(hit[["name", "GeneID", "contig", "l2fc_gam", "padj_gam", "l2fc_col", "padj_col",
           "n_span", "n_focus", "focus_cohorts", "focus_countries", "max_delta_i"]]
      .to_string(index=False) if len(hit) else "  *** NO OBP in the family overlaps any West African H12 signal ***")

print(f"\n  OBPs in family: {len(obp)} | with span overlap: {(obp['n_span']>0).sum()} | "
      f"with focus overlap: {(obp['n_focus']>0).sum()}")
print(f"  OBPs sig-up in BOTH species: {int(obp['sig_both'].sum())}")
print(obp[obp["sig_both"]][["name", "contig", "l2fc_gam", "padj_gam", "l2fc_col", "padj_col",
                            "n_span", "n_focus"]].sort_values("l2fc_gam", ascending=False).to_string(index=False))

print("\n" + "=" * 120)
print("C2. ABC family — focus-level sweep overlap (the strict criterion ABCH1 meets)")
print("=" * 120)
ah = abc[abc["n_focus"] > 0].sort_values("focus_cohorts", ascending=False)
print(ah[["name", "GeneID", "contig", "l2fc_gam", "padj_gam", "l2fc_col", "padj_col",
          "n_span", "n_focus", "focus_cohorts", "focus_countries", "max_delta_i"]].to_string(index=False))

print("\n" + "=" * 120)
print("C3. ABCH1 — does the ABCC 3R array survive the focus-level criterion?")
print("=" * 120)
for gid in ["AGAP002638", "AGAP008436", "AGAP008437", "AGAP027980", "AGAP028128", "AGAP003680"]:
    nm = abc.loc[abc["GeneID"] == gid, "name"]
    print(f"  {gid} ({nm.iloc[0] if len(nm) else '?'}): "
          f"span={len(signal_hits(gid,'span'))}, focus={len(signal_hits(gid,'focus'))}")

# ── Genome-wide context: how exceptional is ABCH1's 8-cohort focus overlap? ──
print("\n" + "=" * 120)
print("C4. Genome-wide context — distribution of focus-level cohort counts per gene")
print("=" * 120)
allf = []
for gid in genes_df["GeneID"]:
    f = signal_hits(gid, "focus")
    allf.append((gid, len(f["cohort_id"].unique()) if len(f) else 0))
fdf = pd.DataFrame(allf, columns=["GeneID", "focus_cohorts"])
print(fdf["focus_cohorts"].value_counts().sort_index().to_string())
n_ge8 = int((fdf["focus_cohorts"] >= 8).sum())
print(f"\n  Genes with focus overlap in >=8 cohorts: {n_ge8} / {len(fdf)} "
      f"({100*n_ge8/len(fdf):.2f}%) -> ABCH1 is in the top {100*n_ge8/len(fdf):.2f}% of the genome")
fdf.merge(genes_df[["GeneID", "contig", "Name"]], on="GeneID").sort_values(
    "focus_cohorts", ascending=False).head(25).to_csv("results/abc_obp/top_swept_genes.csv", index=False)
print("\n  Top 25 most-swept genes genome-wide (focus level):")
print(fdf.merge(genes_df[["GeneID", "contig", "Name"]], on="GeneID")
      .sort_values("focus_cohorts", ascending=False).head(25).to_string(index=False))

obp.to_csv("results/abc_obp/obp_profile.csv", index=False)
abc.to_csv("results/abc_obp/abc_profile.csv", index=False)
