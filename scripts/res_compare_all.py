import sys
from pathlib import Path
import pandas as pd
import yaml

# Identify directories and files
if len(sys.argv) != 2:
    raise ValueError("Usage: python res_compare_all.py <species name directory>")

species = sys.argv[1]
data_directory = Path("data") / species
data_files = sorted(data_directory.glob("*/Resistance.DNA.xlsx"))

# STEP 1: Read resistance thresholds from config.yaml
with open("config.yaml") as file:
    config = yaml.safe_load(file)

min_coverage = config["resistance"]["min_coverage"]
min_depth = config["resistance"]["min_depth"]

# STEP 2: read all variants from the all the samples
samples = dict()
for file in data_files:
    df = pd.read_excel(file)
    sample_name = file.parent.name
    samples[sample_name] = df

# STEP 3: Match resistance genes across samples
variants = dict()
for sample_name, df in samples.items():
    for index, row in df.iterrows():
        key = (row['Gene'], row['Class'], row['Mechanism'], row['Type'])
        if key not in variants:
            variants[key] = dict()
        if sample_name not in variants[key]:
            variants[key][sample_name] = list()
        variants[key][sample_name].append({'Coverage(%)': row['Coverage(%)'], 'Depth': row['Depth']})

# STEP 4: Exclude variants shared by all samples
results = list()
for variant, sample_data in variants.items():
    if len(sample_data) != len(samples):
        results.append({'variant': variant, 'samples': sample_data})

# STEP 5: Apply resistance gene coverage and depth thresholds
filtered_results = list()
for result in results:
    filtered_samples = dict()
    for sample_name, hits in result["samples"].items():
        passing_hits = list()
        for hit in hits:
            coverage = pd.to_numeric(hit["Coverage(%)"], errors="coerce")
            depth = pd.to_numeric(hit["Depth"], errors="coerce")
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

# STEP 6: create data frame
rows = []
for result in results:
    gene, drug_class, mechanism, gene_type = result['variant']
    sample_data = result['samples']
    row = {
        'Gene': gene,
        'Class': drug_class,
        'Mechanism': mechanism,
        'Type': gene_type
    }
    for sample_name in samples:
        hits = sample_data.get(sample_name, [])
        if hits:
            row[f'{sample_name}_Coverage(%)'] = '; '.join(str(hit['Coverage(%)']) for hit in hits)
            row[f'{sample_name}_Depth'] = '; '.join(str(hit['Depth']) for hit in hits)
        else:
            row[f'{sample_name}_Coverage(%)'] = pd.NA
            row[f'{sample_name}_Depth'] = pd.NA
    rows.append(row)
comparison_df = pd.DataFrame(rows)

# STEP 7: Export comparison results
output_directory = Path("results") / species
output_directory.mkdir(parents=True, exist_ok=True)
output_file = output_directory / "res_comparison.xlsx"
comparison_df.to_excel(output_file, index=False)
print(f"Results saved to: {output_file.resolve()}")
