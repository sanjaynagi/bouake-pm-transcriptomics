"""B9 — regenerate the ABC differential-expression and sweep-overlap tables with the strict ABC definition.

The ABC family is the ABC nucleotide-binding domain (PF00005) or an Abc* gene name (see setup.py).
Writes results/abc_obp/bouake_de_abc_{gambiae,coluzzii}.csv (every family gene tested, ranked by
log2 fold change) and results/abc_obp/sweep_overlap_abc_gambiae.csv (the significantly overexpressed
An. gambiae members with their West African H12 sweep overlap), and prints the counts quoted in the
manuscript.
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

from setup import FAM, ann_by_gene, col, gam, genes_df

PFAM = "PF00005"


def family_table(*, de: pd.DataFrame) -> pd.DataFrame:
    """Family members in one DE table, with how each qualified for the family."""
    pfam = ann_by_gene.set_index("gene_id")["pfam"].fillna("")
    table = de[de["GeneID"].isin(FAM["ABC"])].copy()
    has_domain = table["GeneID"].map(lambda g: PFAM in pfam.get(g, ""))
    has_name = table["GeneName"].fillna("").str.contains(r"^Abc", case=False, regex=True)
    table["family"] = "ABC"
    table["family_source"] = np.select([has_domain & has_name, has_domain], ["pfam+name", "pfam"], default="name")
    table["pfam_domains"] = table["GeneID"].map(pfam)
    return table.sort_values("log2FoldChange", ascending=False)


def significant_counts(*, table: pd.DataFrame) -> tuple[int, int]:
    significant = table[table["padj"] < 0.05]
    return int((significant["log2FoldChange"] > 0).sum()), int((significant["log2FoldChange"] < 0).sum())


if __name__ == "__main__":
    tables = {"gambiae": family_table(de=gam), "coluzzii": family_table(de=col)}
    for species, table in tables.items():
        table.to_csv(f"results/abc_obp/bouake_de_abc_{species}.csv", index=False)
    up_gam, down_gam = significant_counts(table=tables["gambiae"])
    up_col, down_col = significant_counts(table=tables["coluzzii"])
    print(f"family size {len(FAM['ABC'])}; tested in gambiae {len(tables['gambiae'])}, coluzzii {len(tables['coluzzii'])}")
    print(f"gambiae: {up_gam} up, {down_gam} down; coluzzii: {up_col} up, {down_col} down")
    print("Fisher exact, up vs down, gambiae vs coluzzii:",
          fisher_exact([[up_gam, down_gam], [up_col, down_col]]))

    sweeps = (tables["gambiae"][(tables["gambiae"]["padj"] < 0.05) & (tables["gambiae"]["log2FoldChange"] > 0)]
              .merge(genes_df[["GeneID", "contig", "start", "end", "n_sweeps", "max_delta_i", "countries", "cohorts"]],
                     on="GeneID", how="left"))
    sweeps.to_csv("results/abc_obp/sweep_overlap_abc_gambiae.csv", index=False)
