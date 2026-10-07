import pandas as pd

# Input and output file names
input_file = "BouakePM_diffexp.xlsx"
output_subset_file = "BouakePM_diffexp_subset.xlsx"
output_families_file = "BouakePM_diffexp_P450_ABC_GST_COE.xlsx"

# List of target genes
target_genes = [
    "AGAP001356", "AGAP006227", "AGAP006228", "AGAP006723", "AGAP006724",
    "AGAP006725", "AGAP006726", "AGAP006727", "AGAP006728", "AGAP008212",
    "AGAP008213", "AGAP008214", "AGAP005373", "AGAP005372", "AGAP008052",
    "AGAP000818",
]

# Keywords for second filtering pass
search_terms = ["cytochrome P450", "ABC transporter", "glutathione S", "esterase"]

columns_to_drop = ["Log2FoldChange", "lfcSE", "stat", "pvalue"]

# Read all sheets
xlsx = pd.ExcelFile(input_file)

# Create a writer for each output Excel file
with pd.ExcelWriter(output_subset_file, engine="openpyxl") as writer_subset, \
     pd.ExcelWriter(output_families_file, engine="openpyxl") as writer_families:

    for sheet_name in xlsx.sheet_names:
        print(f"Processing sheet: {sheet_name}")
        df = pd.read_excel(xlsx, sheet_name=sheet_name)

        # Normalize column names (important in case they vary by sheet)
        df.columns = df.columns.str.strip()

        # Filter 1: Specific genes
        df_subset = df[df["GeneID"].isin(target_genes)]

        # Filter 2: GeneDescription keywords (case-insensitive)
        pattern = "|".join(search_terms)
        df_families = df[
            df["GeneDescription"]
            .astype(str)
            .str.contains(pattern, case=False, na=False)
        ]

        df_subset = df_subset.sort_values(by='absolute_diff', ascending=False)
        df_families = df_families.sort_values(by='absolute_diff', ascending=False)

        # Write each filtered dataframe to its respective Excel
        df_subset.drop(columns=columns_to_drop, errors='ignore').to_excel(writer_subset, sheet_name=sheet_name, index=False)
        df_families.to_excel(writer_families, sheet_name=sheet_name, index=False)

print("✅ Done! Subset and family-filtered Excel files created.")