"""Differential-expression values quoted in the Results, taken from the RNA-Seq-Pop DESeq2 tables.

Reads results/BouakePM_diffexp.xlsx (one sheet per contrast) and writes
results/rnaseq/de_numbers.csv with the fold change (2 ** log2FoldChange) and adjusted P of each gene
named in the text, so each number in the manuscript can be checked against its source.
"""

import pandas as pd

from bouake.paths import RESULTS

SHEETS = {"gambiae": "gambiaeCont_gambiaePM", "coluzzii": "coluzziiCont_coluzziiPM"}
GENES = {  # label -> gene name or id, in the order quoted in the text
    "gambiae": ["CYP6M2", "CYP6M3", "CYP6M4", "COEAE3G", "COEAE7G", "COEBE1C", "COEBE3C", "SAP2",
                "CYP9K1", "Ace1", "ABCH1"],
    "coluzzii": ["ABCH1"],
}


def quoted_values(*, species: str) -> pd.DataFrame:
    """Fold change and adjusted P for the genes named in the text, for one species' survivor contrast."""
    table = pd.read_excel(RESULTS / "BouakePM_diffexp.xlsx", sheet_name=SHEETS[species])
    name = table["GeneName"].fillna("").str.upper()
    rows = []
    for gene in GENES[species]:
        hit = table[(name == gene.upper()) | (table["GeneID"] == gene)]
        assert len(hit) == 1, (species, gene, len(hit))
        rows.append(dict(species=species, gene=gene, fold_change=2 ** hit["log2FoldChange"].iloc[0],
                         padj=hit["padj"].iloc[0]))
    return pd.DataFrame(rows)


if __name__ == "__main__":
    (RESULTS / "rnaseq").mkdir(exist_ok=True)
    values = pd.concat([quoted_values(species=s) for s in SHEETS])
    values.to_csv(RESULTS / "rnaseq" / "de_numbers.csv", index=False)
    print(values.round(4).to_string(index=False))
