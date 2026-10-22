from pathlib import Path
import sys
import pandas as pd
import yaml

# Species supplied by Snakemake
if len(sys.argv) != 2:
    raise ValueError("Usage: python cnv_compare_all.py <species>")

species = sys.argv[1]
data_dir = Path("data") / species

with open("config.yaml") as file:
    config = yaml.safe_load(file)

cnv_config = config["cnv"]
threshold = cnv_config["min_reciprocal_overlap"]


def is_match(a, b):
    if a["Chromosome"] != b["Chromosome"]:
        return False

    if a["Type"] != b["Type"]:
        return False

    overlap = max(
        0,
        min(a["End"], b["End"]) - max(a["Start"], b["Start"]) + 1
    )

    length_a = a["End"] - a["Start"] + 1
    length_b = b["End"] - b["Start"] + 1

    return (
        overlap / length_a >= threshold
        and overlap / length_b >= threshold
    )


# Discover samples and read all CNV calls before filtering
paths = sorted(data_dir.glob("*/Anno_cnv.xls"))

if not paths:
    raise FileNotFoundError("No CNV files found.")

sample_names = [path.parent.name for path in paths]
calls = []

for path in paths:
    df = pd.read_csv(path, sep="\t")

    # Keep all original CNV calls for cross-sample comparison
    for _, row in df.iterrows():
        calls.append({
            "Sample": path.parent.name,
            "Chromosome": row["#Chr"],
            "Start": int(row["Pos1"]),
            "End": int(row["Pos2"]),
            "Type": row["Type"],
            "Normalized_RD": row["Normalized_RD"]
        })


# Group overlapping CNVs
calls.sort(
    key=lambda x: (
        x["Chromosome"],
        x["Type"],
        x["Start"],
        x["End"]
    )
)

groups = []

for call in calls:
    matching_groups = [
        group
        for group in groups
        if all(is_match(call, member) for member in group)
    ]

    if len(matching_groups) == 1:
        matching_groups[0].append(call)

    elif len(matching_groups) == 0:
        groups.append([call])

    else:
        raise ValueError(
            "Ambiguous CNV grouping: one call matches multiple groups. "
            "Manual review is required."
        )


# Retain groups absent from at least one sample
results = []

for group in groups:
    present_in = sorted({call["Sample"] for call in group})

    # Exclude CNVs reported in all four samples
    if len(present_in) == len(sample_names):
        continue

    # Apply RD thresholds only after cross-sample comparison
    passing_calls = [
        call for call in group
        if (
            call["Type"] == "deletion"
            and call["Normalized_RD"] < cnv_config["deletion_max_rd"]
        )
        or (
            call["Type"] == "duplication"
            and call["Normalized_RD"] > cnv_config["duplication_min_rd"]
        )
    ]

    # Skip CNV groups with no calls passing the RD thresholds
    if not passing_calls:
        continue

    rd_by_sample = "; ".join(
        f'{call["Sample"]}: {call["Normalized_RD"]:.6f}'
        for call in sorted(group, key=lambda x: x["Sample"])
    )

    results.append({
        "Chromosome": group[0]["Chromosome"],
        "Start": min(call["Start"] for call in group),
        "End": max(call["End"] for call in group),
        "Type": group[0]["Type"],
        "Present_in_samples": ", ".join(present_in),
        "Normalized_RD_by_sample": rd_by_sample
    })


columns = [
    "Chromosome",
    "Start",
    "End",
    "Type",
    "Present_in_samples",
    "Normalized_RD_by_sample"
]

result_df = pd.DataFrame(results, columns=columns)

result_df = result_df.sort_values(
    ["Chromosome", "Start", "End"]
).reset_index(drop=True)

print("Distinct CNV events requiring review:", len(result_df))

if not result_df.empty:
    print(result_df.to_string(index=False))

output = Path("results") / species / "cnv_comparison_all.csv"

output.parent.mkdir(parents=True, exist_ok=True)

result_df.to_csv(output, index=False)

print(f"\nSaved to {output}")


