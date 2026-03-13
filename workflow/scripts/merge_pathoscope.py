"""
Merge per-sample Pathoscope SAM reports into a combined TSV with taxonomy resolution.
Called by Snakemake; uses snakemake.input and snakemake.output.
"""
import pandas as pd
import os

try:
    from ete3 import NCBITaxa
    ncbi = NCBITaxa()
    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False
    print("WARNING: ete3 not available — taxonomy columns will not be resolved.")

TAXA_RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
BASE_COLS  = ["tid", "superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "name"]

input_files = list(snakemake.input)
output_file = snakemake.output[0]

combined = pd.DataFrame(columns=BASE_COLS)

for fpath in input_files:
    sample_id = os.path.basename(fpath).split("-sam-report")[0]
    try:
        with open(fpath) as f:
            first_line = f.readline()
        df = pd.read_csv(fpath, sep="\t", skiprows=[0])
        df = df[df["Genome"].str.contains("ti", na=False)]
    except Exception as e:
        print(f"WARNING: Could not read {fpath}: {e}")
        continue

    rows = []
    for _, row in df.iterrows():
        parts = str(row["Genome"]).split("|")
        if parts[0] != "ti":
            continue
        tid = parts[1]
        record = {col: "Unknown" for col in BASE_COLS}
        record["tid"] = tid
        record[sample_id] = row["Final Guess"]

        if HAS_ETE3:
            try:
                lineage = ncbi.get_lineage(tid)
                translator = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(list(translator.keys()))
                for taxid, rank in ranks.items():
                    if rank in TAXA_RANKS:
                        record[rank] = translator[taxid]
                record["name"] = "; ".join(
                    record.get(r, "Unknown") for r in TAXA_RANKS
                )
            except Exception:
                pass

        rows.append(record)

    if rows:
        sample_df = pd.DataFrame(rows)
        combined = pd.concat([combined, sample_df], axis=0, ignore_index=True)

# Deduplicate by tid, merge abundance columns
if not combined.empty:
    combined = combined.groupby("tid").first().reset_index()

os.makedirs(os.path.dirname(output_file), exist_ok=True)
combined.to_csv(output_file, sep="\t", index=False)
print(f"Pathoscope merge complete: {len(combined)} taxa written to {output_file}")
