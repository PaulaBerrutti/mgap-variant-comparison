from pathlib import Path
import sys
import pandas as pd
import yaml

# Species supplied by Snakemake
if len(sys.argv) != 2:
    raise ValueError("Usage: python sv_compare_all.py <species>")

species = sys.argv[1]
data_dir = Path("data") / species

# Read filtering thresholds
with open("config.yaml") as file:
    config = yaml.safe_load(file)

sv_config = config["sv"]

# Discover SV files
paths = sorted(data_dir.glob("*/Anno_sv.xls"))

if not paths:
    raise FileNotFoundError("No MGAP SV files found.")

samples = {}
variants = {}

required = [
    "#Chr1", "Pos1", "Chr2", "Pos2",
    "Type", "Num_Reads", "Score",
    "Estimated_allele_frequency"
]

# Read and filter each sample
for path in paths:

    df = pd.read_csv(path, sep="\t")

    missing = [col for col in required if col not in df.columns]

    if missing:
        raise ValueError(f"{path}: missing columns {missing}")

    sample_name = path.parent.name

    # Keep all calls for comparison before filtering
    filtered = df.copy()

    # Register the sample for cross-sample comparison
    samples[sample_name] = filtered

    # Collect matching SVs
    for _, row in filtered.iterrows():

        key = (
            row["#Chr1"],
            int(row["Pos1"]),
            row["Chr2"],
            int(row["Pos2"]),
            row["Type"]
        )

        if key not in variants:
            variants[key] = {
                "samples": set(),
                "evidence": []
            }

        variants[key]["samples"].add(sample_name)

        variants[key]["evidence"].append({
            "Sample": sample_name,
            "Num_Reads": row["Num_Reads"],
            "Score": row["Score"],
            "Estimated_allele_frequency":
                row["Estimated_allele_frequency"]
        })


# Keep SVs not present in all samples
results = []

for key, data in variants.items():

    present_in = data["samples"]

    if len(present_in) < len(samples):

        chr1, pos1, chr2, pos2, sv_type = key

        for evidence in data["evidence"]:

            # Apply score threshold only after comparison
            if evidence["Score"] < sv_config["min_score"]:
                continue


            results.append({
                "Chromosome_1": chr1,
                "Position_1": pos1,
                "Chromosome_2": chr2,
                "Position_2": pos2,
                "Type": sv_type,
                "Num_Reads": evidence["Num_Reads"],
                "Score": evidence["Score"],
                "Estimated_allele_frequency":
                    evidence["Estimated_allele_frequency"],
                "Present_in_samples": ", ".join(sorted(present_in))
            })


columns = [
    "Chromosome_1",
    "Position_1",
    "Chromosome_2",
    "Position_2",
    "Type",
    "Num_Reads",
    "Score",
    "Estimated_allele_frequency",
    "Present_in_samples"
]

result_df = pd.DataFrame(results, columns=columns)

result_df = result_df.sort_values(
    ["Chromosome_1", "Position_1"]
).reset_index(drop=True)

print("SV calls requiring review:", len(result_df))

if not result_df.empty:
    print(result_df.to_string(index=False))

# Save results
output = Path("results") / species / "sv_comparison_all.csv"

output.parent.mkdir(parents=True, exist_ok=True)

result_df.to_csv(output, index=False)

print(f"\nSaved to {output}")
