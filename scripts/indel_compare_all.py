from pathlib import Path
import sys
import pandas as pd

# Species supplied by Snakemake
if len(sys.argv) != 2:
    raise ValueError("Usage: python indel_compare_all.py <species>")

species = sys.argv[1]
data_dir = Path("data") / species

# Find indel files
paths = sorted(data_dir.glob("*/variant.indel.ann.tab.xlsx"))

if not paths:
    raise FileNotFoundError("No MGAP indel files found.")

samples = {}

# Read indel records from every sample
for path in paths:

    df = pd.read_excel(path)

    required = ["CHROM", "POS", "TYPE", "REF", "MAF"]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    # Keep all records for comparison before filtering
    samples[path.parent.name] = df


# Group identical indel calls across samples
variants = {}

for sample_name, df in samples.items():

    for _, row in df.iterrows():

        key = (
            row["CHROM"],
            row["POS"],
            row["TYPE"],
            row["REF"],
            row["MAF"]
        )

        if key not in variants:
            variants[key] = set()

        variants[key].add(sample_name)


# Keep indels not reported in all samples
results = []



for key, present_in in variants.items():

    # Exclude variants reported in all four samples
    if len(present_in) == len(samples):
        continue

    chrom, pos, variant_type, ref, maf = key

    # Apply variant-type selection after comparison
    if variant_type not in ["ins", "del"]:
        continue

    results.append({
            "Chromosome": chrom,
            "Position": pos,
            "Type": variant_type,
            "REF": ref,
            "MAF": maf,
            "Present_in_samples": ", ".join(sorted(present_in))
        })


# Create results table
columns = [
    "Chromosome",
    "Position",
    "Type",
    "REF",
    "MAF",
    "Present_in_samples"
]

result_df = pd.DataFrame(results, columns=columns)

result_df = result_df.sort_values(
    ["Chromosome", "Position"]
).reset_index(drop=True)


# Display results
print("Indel variants requiring review:", len(result_df))

if not result_df.empty:
    print(result_df.to_string(index=False))


# Save results
output = Path("results") / species / "indel_comparison_all.csv"

output.parent.mkdir(parents=True, exist_ok=True)

result_df.to_csv(output, index=False)

print(f"\nSaved to {output}")
