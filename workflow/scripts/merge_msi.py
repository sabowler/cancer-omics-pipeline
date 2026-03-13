"""
Merge per-sample MSI scores into a single Excel report.
Called by Snakemake; uses snakemake.input and snakemake.output.
"""
import pandas as pd
import os

input_files = snakemake.input
output_file = snakemake.output[0]

records = []
for fpath in input_files:
    sample = os.path.basename(fpath)
    try:
        df = pd.read_csv(fpath, sep="\t")
        records.append({
            "Sample":                  sample,
            "Total_Number_of_Sites":   df["Total_Number_of_Sites"][0],
            "Number_of_Somatic_Sites": df["Number_of_Somatic_Sites"][0],
            "MSI_Score":               df["%"][0],
        })
    except Exception as e:
        print(f"WARNING: Could not parse {fpath}: {e}")

result = pd.DataFrame(records)
os.makedirs(os.path.dirname(output_file), exist_ok=True)
result.to_excel(output_file, index=False)
print(f"MSI merge complete: {len(result)} samples written to {output_file}")
