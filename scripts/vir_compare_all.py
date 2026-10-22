
import sys
from pathlib import Path
import pandas as pd
import yaml


# ==========================================
# Identify species and input files
# ==========================================

if len(sys.argv) != 2:
    raise ValueError("Usage: python vir_compare_all.py <species>")

species = sys.argv[1]

data_directory = Path("data") / species

data_files = sorted(
    data_directory.glob("*/Virulence.DNA.xlsx")
)

if not data_files:
    raise FileNotFoundError(
        f"No Virulence.DNA.xlsx files found for {species}"
    )


# ==========================================
# STEP 1: Read thresholds from config.yaml
# ==========================================

with open("config.yaml") as file:
    config = yaml.safe_load(file)

min_coverage = config["virulence"]["min_coverage"]
min_depth = config["virulence"]["min_depth"]


# ==========================================
# STEP 2: Read virulence genes from all samples
# ==========================================

samples = {}

required_columns = [
    "Virulence Factor",
    "Pathogen",
    "VF Protein",
    "Coverage(%)",
    "Depth"
]

for file in data_files:

    df = pd.read_excel(file)

    missing = [
        col for col in required_columns
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{file}: missing columns {missing}"
        )

    sample_name = file.parent.name

    samples[sample_name] = df


# ==========================================
# STEP 3: Match virulence genes across samples
# ==========================================

variants = {}

for sample_name, df in samples.items():

    for _, row in df.iterrows():

        key = (
            row["Virulence Factor"],
            row["VF Protein"]
        )

        if key not in variants:
            variants[key] = {}

        if sample_name not in variants[key]:
            variants[key][sample_name] = []

        variants[key][sample_name].append({
            "Pathogen": row["Pathogen"],
            "Coverage(%)": row["Coverage(%)"],
            "Depth": row["Depth"]
        })


# ==========================================
# STEP 4: Exclude findings shared by all samples
# ==========================================

results = []

for variant, sample_data in variants.items():

    if len(sample_data) != len(samples):

        results.append({
            "variant": variant,
            "samples": sample_data
        })


# ==========================================
# STEP 5: Apply coverage and depth thresholds
# ==========================================

filtered_results = []

for result in results:

    filtered_samples = {}

    for sample_name, hits in result["samples"].items():

        passing_hits = []

        for hit in hits:

            coverage = pd.to_numeric(
                str(hit["Coverage(%)"]).replace("%", ""),
                errors="coerce"
            )

            depth = pd.to_numeric(
                hit["Depth"],
                errors="coerce"
            )

            if coverage >= min_coverage and depth >= min_depth:
                passing_hits.append(hit)

        if passing_hits:
            filtered_samples[sample_name] = passing_hits

    if filtered_samples:

        filtered_results.append({
            "variant": result["variant"],
            "samples": filtered_samples
        })

results = filtered_results


# ==========================================
# STEP 6: Create comparison DataFrame
# ==========================================

rows = []

for result in results:

    virulence_factor, vf_protein = result["variant"]

    sample_data = result["samples"]

    row = {
        "Virulence Factor": virulence_factor,
        "VF Protein": vf_protein
    }

    for sample_name in samples:

        hits = sample_data.get(sample_name, [])

        if hits:

            row[f"{sample_name}_Pathogen"] = "; ".join(
                str(hit["Pathogen"]) for hit in hits
            )

            row[f"{sample_name}_Coverage(%)"] = "; ".join(
                str(hit["Coverage(%)"]) for hit in hits
            )

            row[f"{sample_name}_Depth"] = "; ".join(
                str(hit["Depth"]) for hit in hits
            )

        else:

            row[f"{sample_name}_Pathogen"] = pd.NA
            row[f"{sample_name}_Coverage(%)"] = pd.NA
            row[f"{sample_name}_Depth"] = pd.NA

    rows.append(row)


# Keep column headers even when no findings pass filtering
columns = ["Virulence Factor", "VF Protein"]

for sample_name in samples:
    columns.extend([
        f"{sample_name}_Pathogen",
        f"{sample_name}_Coverage(%)",
        f"{sample_name}_Depth"
    ])

comparison_df = pd.DataFrame(rows, columns=columns)


# ==========================================
# STEP 7: Export comparison results
# ==========================================

output_directory = Path("results") / species

output_directory.mkdir(
    parents=True,
    exist_ok=True
)

output_file = output_directory / "vir_comparison.xlsx"

comparison_df.to_excel(
    output_file,
    index=False
)

print(
    f"Virulence findings requiring review: {len(comparison_df)}"
)

print(f"Results saved to: {output_file.resolve()}")
