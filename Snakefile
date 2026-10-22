
# ==========================================
# MGAP Variant Comparison Pipeline
# Multiple bacterial species
# ==========================================

configfile: "config.yaml"

from pathlib import Path

data_dir = Path("data")


# ==========================================
# Files required for each sample
# ==========================================

required_files = [
    "Anno_cnv.xls",
    "Anno_sv.xls",
    "variant.snp.ann.tab.xlsx",
    "variant.indel.ann.tab.xlsx",
    "Resistance.DNA.xlsx",
    "Virulence.DNA.xlsx"
]


# ==========================================
# Discover species and samples
# ==========================================

species_samples = {}

for species_dir in sorted(data_dir.iterdir()):

    if not species_dir.is_dir():
        continue

    sample_dirs = sorted(
        folder for folder in species_dir.iterdir()
        if folder.is_dir()
    )

    # Skip species without sample folders
    if not sample_dirs:
        print(
            f"Skipping {species_dir.name}: no samples yet."
        )
        continue

    # Check required MGAP files
    for sample_dir in sample_dirs:

        for filename in required_files:

            filepath = sample_dir / filename

            if not filepath.is_file():
                raise FileNotFoundError(
                    f"Missing required file: {filepath}"
                )

    species_samples[species_dir.name] = [
        folder.name for folder in sample_dirs
    ]


if not species_samples:
    raise ValueError("No species with samples found.")


species_names = sorted(species_samples)

for species, samples in species_samples.items():

    print(
        f"Found {species}: {len(samples)} samples "
        f"({', '.join(samples)})"
    )


# ==========================================
# Helper: files belonging to one species
# ==========================================

def sample_files(wildcards, filename):

    return [
        str(
            data_dir
            / wildcards.species
            / sample
            / filename
        )
        for sample in species_samples[wildcards.species]
    ]


# ==========================================
# Final output files
# ==========================================

reports = [
    "variant_filter_summary.csv",
    "cnv_comparison_all.csv",
    "snp_comparison_all.csv",
    "indel_comparison_all.csv",
    "sv_comparison_all.csv",
    "res_comparison.xlsx",
    "vir_comparison.xlsx"
]


rule all:
    input:
        expand(
            "results/{species}/{report}",
            species=species_names,
            report=reports
        )


# ==========================================
# Variant filtering summary
# ==========================================

rule variant_filter_summary:
    input:
        cnv_files=lambda wc: sample_files(
            wc, "Anno_cnv.xls"
        ),
        sv_files=lambda wc: sample_files(
            wc, "Anno_sv.xls"
        ),
        snp_files=lambda wc: sample_files(
            wc, "variant.snp.ann.tab.xlsx"
        ),
        indel_files=lambda wc: sample_files(
            wc, "variant.indel.ann.tab.xlsx"
        ),
        res_files=lambda wc: sample_files(
            wc, "Resistance.DNA.xlsx"
        ),
        vir_files=lambda wc: sample_files(
            wc, "Virulence.DNA.xlsx"
        ),
        config_file="config.yaml",
        script_file="scripts/variant_filter_summary.py",
        res_script="scripts/res_compare_all.py",
        vir_script="scripts/vir_compare_all.py"
    output:
        "results/{species}/variant_filter_summary.csv"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# CNV comparison
# ==========================================

rule cnv_compare:
    input:
        cnv_files=lambda wc: sample_files(
            wc, "Anno_cnv.xls"
        ),
        config_file="config.yaml",
        script_file="scripts/cnv_compare_all.py"
    output:
        "results/{species}/cnv_comparison_all.csv"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# SNP/MNP/complex comparison
# ==========================================

rule snp_compare:
    input:
        snp_files=lambda wc: sample_files(
            wc, "variant.snp.ann.tab.xlsx"
        ),
        script_file="scripts/snp_compare_all.py"
    output:
        "results/{species}/snp_comparison_all.csv"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# Indel comparison
# ==========================================

rule indel_compare:
    input:
        indel_files=lambda wc: sample_files(
            wc, "variant.indel.ann.tab.xlsx"
        ),
        script_file="scripts/indel_compare_all.py"
    output:
        "results/{species}/indel_comparison_all.csv"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# SV comparison
# ==========================================

rule sv_compare:
    input:
        sv_files=lambda wc: sample_files(
            wc, "Anno_sv.xls"
        ),
        config_file="config.yaml",
        script_file="scripts/sv_compare_all.py"
    output:
        "results/{species}/sv_comparison_all.csv"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# Resistance gene comparison
# ==========================================

rule res_compare:
    input:
        res_files=lambda wc: sample_files(
            wc, "Resistance.DNA.xlsx"
        ),
        config_file="config.yaml",
        script_file="scripts/res_compare_all.py"
    output:
        "results/{species}/res_comparison.xlsx"
    shell:
        "python {input.script_file} {wildcards.species}"


# ==========================================
# Virulence gene comparison
# ==========================================

rule vir_compare:
    input:
        vir_files=lambda wc: sample_files(
            wc, "Virulence.DNA.xlsx"
        ),
        config_file="config.yaml",
        script_file="scripts/vir_compare_all.py"
    output:
        "results/{species}/vir_comparison.xlsx"
    shell:
        "python {input.script_file} {wildcards.species}"
