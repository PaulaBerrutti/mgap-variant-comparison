from pathlib import Path
import sys
import pandas as pd

# Identify species and data
if len(sys.argv) != 2:
    raise ValueError("Usage: python snp_compare_all.py <species>")

species = sys.argv[1]
data_dir = Path("data") / species

# Discover MGAP SNP files
paths = sorted(data_dir.glob("*/variant.snp.ann.tab.xlsx"))

if not paths:
    raise FileNotFoundError("No MGAP SNP files found.")

# STEP 1: Read all variant calls from all samples
samples = {}

for path in paths:
    df = pd.read_excel(path)

    required = ["CHROM", "POS", "TYPE", "REF", "MAF"]
    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    sample_name = path.parent.name
    samples[sample_name] = df

# STEP 2: Match variants across samples
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

# STEP 3: Exclude variants shared by all samples
# STEP 4: Select SNP, MNP and complex variants
results = []

for key, present_in in variants.items():

    if len(present_in) == len(samples):
        continue

    chrom, pos, variant_type, ref, maf = key

    if variant_type not in ["snp", "mnp", "complex"]:
        continue

    # STEP 5: Retain variants for manual investigation
    results.append({
        "Chromosome": chrom,
        "Position": pos,
        "Type": variant_type,
        "REF": ref,
        "MAF": maf,
        "Present_in_samples": ", ".join(sorted(present_in))
    })

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

print("SNP/MNP/complex variants requiring review:", len(result_df))

if not result_df.empty:
    print(result_df.to_string(index=False))

# Save results
output = Path("results") / species / "snp_comparison_all.csv"

output.parent.mkdir(parents=True, exist_ok=True)

result_df.to_csv(output, index=False)

print(f"\nSaved to {output}")


