
from pathlib import Path
import sys
import pandas as pd
import yaml


# ==========================================
# Identify species and configuration
# ==========================================

if len(sys.argv) != 2:
    raise ValueError(
        "Usage: python variant_filter_summary.py <species>"
    )

species = sys.argv[1]

with open("config.yaml") as file:
    config = yaml.safe_load(file)

data_dir = Path("data") / species

output = (
    Path("results")
    / species
    / "variant_filter_summary.csv"
)


# ==========================================
# Discover samples
# ==========================================

samples = sorted(
    folder
    for folder in data_dir.iterdir()
    if folder.is_dir()
    and (folder / "Anno_cnv.xls").is_file()
)

if not samples:
    raise FileNotFoundError(
        "No sample directories found."
    )

sample_names = {folder.name for folder in samples}

results = []


# ==========================================
# Helper: add a summary row
# ==========================================

def add_summary(
    variant_type,
    original,
    non_shared,
    final_calls,
    final_events,
    method
):
    results.append({
        "Variant_type": variant_type,
        "Original_calls": original,
        "Non_shared_calls": non_shared,
        "Removed_by_filter": non_shared - final_calls,
        "Final_calls": final_calls,
        "Final_events": final_events,
        "Filtering_method": method
    })


# ==========================================
# SNP and Indel comparisons
# ==========================================

def summarize_exact_variants(
    filename,
    allowed_types,
    label
):
    variants = {}
    original = 0

    for sample_dir in samples:

        df = pd.read_excel(sample_dir / filename)

        original += len(df)

        for _, row in df.iterrows():

            key = (
                row["CHROM"],
                row["POS"],
                row["TYPE"],
                row["REF"],
                row["MAF"]
            )

            variants.setdefault(
                key, set()
            ).add(sample_dir.name)

    # Exclude variants shared by all samples
    non_shared = {
        key: present
        for key, present in variants.items()
        if present != sample_names
    }

    non_shared_calls = sum(
        len(present)
        for present in non_shared.values()
    )

    # Select allowed variant types
    retained = {
        key: present
        for key, present in non_shared.items()
        if key[2] in allowed_types
    }

    final_calls = sum(
        len(present)
        for present in retained.values()
    )

    add_summary(
        label,
        original,
        non_shared_calls,
        final_calls,
        len(retained),
        f"TYPE in {sorted(allowed_types)}"
    )


summarize_exact_variants(
    "variant.snp.ann.tab.xlsx",
    {"snp", "mnp", "complex"},
    "SNP"
)

summarize_exact_variants(
    "variant.indel.ann.tab.xlsx",
    {"ins", "del"},
    "Indel"
)


# ==========================================
# SV comparison
# ==========================================

sv_variants = {}
sv_original = 0

for sample_dir in samples:

    df = pd.read_csv(
        sample_dir / "Anno_sv.xls",
        sep="\t"
    )

    sv_original += len(df)

    for _, row in df.iterrows():

        key = (
            row["#Chr1"],
            int(row["Pos1"]),
            row["Chr2"],
            int(row["Pos2"]),
            row["Type"]
        )

        sv_variants.setdefault(key, []).append({
            "Sample": sample_dir.name,
            "Score": row["Score"]
        })


sv_non_shared = {
    key: calls
    for key, calls in sv_variants.items()
    if {call["Sample"] for call in calls} != sample_names
}

sv_non_shared_calls = sum(
    len(calls)
    for calls in sv_non_shared.values()
)


sv_retained = {
    key: [
        call
        for call in calls
        if call["Score"] >= config["sv"]["min_score"]
    ]
    for key, calls in sv_non_shared.items()
}

sv_retained = {
    key: calls
    for key, calls in sv_retained.items()
    if calls
}


add_summary(
    "SV",
    sv_original,
    sv_non_shared_calls,
    sum(
        len(calls)
        for calls in sv_retained.values()
    ),
    len(sv_retained),
    f"BreakDancer Score >= {config['sv']['min_score']}"
)


# ==========================================
# CNV comparison
# ==========================================

cnv_calls = []
cnv_original = 0

overlap_threshold = (
    config["cnv"]["min_reciprocal_overlap"]
)


def is_match(a, b):

    if a["Chromosome"] != b["Chromosome"]:
        return False

    if a["Type"] != b["Type"]:
        return False

    overlap = max(
        0,
        min(a["End"], b["End"])
        - max(a["Start"], b["Start"])
        + 1
    )

    length_a = a["End"] - a["Start"] + 1
    length_b = b["End"] - b["Start"] + 1

    return (
        overlap / length_a >= overlap_threshold
        and overlap / length_b >= overlap_threshold
    )


for sample_dir in samples:

    df = pd.read_csv(
        sample_dir / "Anno_cnv.xls",
        sep="\t"
    )

    cnv_original += len(df)

    for _, row in df.iterrows():

        cnv_calls.append({
            "Sample": sample_dir.name,
            "Chromosome": row["#Chr"],
            "Start": int(row["Pos1"]),
            "End": int(row["Pos2"]),
            "Type": row["Type"],
            "Normalized_RD": row["Normalized_RD"]
        })


cnv_calls.sort(
    key=lambda x: (
        x["Chromosome"],
        x["Type"],
        x["Start"],
        x["End"]
    )
)


cnv_groups = []

for call in cnv_calls:

    matching_groups = [
        group
        for group in cnv_groups
        if all(
            is_match(call, member)
            for member in group
        )
    ]

    if len(matching_groups) == 1:

        matching_groups[0].append(call)

    elif len(matching_groups) == 0:

        cnv_groups.append([call])

    else:

        raise ValueError(
            "Ambiguous CNV grouping: manual review required."
        )


cnv_non_shared = [
    group
    for group in cnv_groups
    if {
        call["Sample"]
        for call in group
    } != sample_names
]


cnv_non_shared_calls = sum(
    len(group)
    for group in cnv_non_shared
)


cnv_retained = []

for group in cnv_non_shared:

    passing_calls = [
        call
        for call in group
        if (
            call["Type"] == "deletion"
            and call["Normalized_RD"]
            < config["cnv"]["deletion_max_rd"]
        )
        or (
            call["Type"] == "duplication"
            and call["Normalized_RD"]
            > config["cnv"]["duplication_min_rd"]
        )
    ]

    if passing_calls:
        cnv_retained.append(passing_calls)


add_summary(
    "CNV",
    cnv_original,
    cnv_non_shared_calls,
    sum(
        len(group)
        for group in cnv_retained
    ),
    len(cnv_retained),
    "Normalized RD thresholds after reciprocal-overlap matching"
)


# ==========================================
# Resistance gene comparison
# ==========================================

res_variants = {}
res_original = 0

for sample_dir in samples:

    df = pd.read_excel(
        sample_dir / "Resistance.DNA.xlsx"
    )

    res_original += len(df)

    for _, row in df.iterrows():

        key = (
            row["Gene"],
            row["Class"],
            row["Mechanism"],
            row["Type"]
        )

        if key not in res_variants:
            res_variants[key] = {}

        if sample_dir.name not in res_variants[key]:
            res_variants[key][sample_dir.name] = []

        res_variants[key][sample_dir.name].append({
            "Coverage(%)": row["Coverage(%)"],
            "Depth": row["Depth"]
        })


# Exclude genes shared by all samples
res_non_shared = {
    key: sample_data
    for key, sample_data in res_variants.items()
    if len(sample_data) != len(samples)
}


res_non_shared_calls = sum(
    len(hits)
    for sample_data in res_non_shared.values()
    for hits in sample_data.values()
)


# Apply resistance thresholds
res_retained = {}

for key, sample_data in res_non_shared.items():

    filtered_samples = {}

    for sample_name, hits in sample_data.items():

        passing_hits = []

        for hit in hits:

            coverage = pd.to_numeric(
                hit["Coverage(%)"],
                errors="coerce"
            )

            depth = pd.to_numeric(
                hit["Depth"],
                errors="coerce"
            )

            if (
                coverage >= config["resistance"]["min_coverage"]
                and depth >= config["resistance"]["min_depth"]
            ):
                passing_hits.append(hit)

        if passing_hits:
            filtered_samples[sample_name] = passing_hits

    if filtered_samples:
        res_retained[key] = filtered_samples


res_final_calls = sum(
    len(hits)
    for sample_data in res_retained.values()
    for hits in sample_data.values()
)


add_summary(
    "Resistance",
    res_original,
    res_non_shared_calls,
    res_final_calls,
    len(res_retained),
    (
        f"Coverage >= {config['resistance']['min_coverage']}% "
        f"and Depth >= {config['resistance']['min_depth']}"
    )
)


# ==========================================
# Virulence gene comparison
# ==========================================

vir_variants = {}
vir_original = 0


# Read all virulence hits
for sample_dir in samples:

    df = pd.read_excel(
        sample_dir / "Virulence.DNA.xlsx"
    )

    vir_original += len(df)

    required_columns = [
        "Virulence Factor",
        "Pathogen",
        "VF Protein",
        "Coverage(%)",
        "Depth"
    ]

    missing = [
        col
        for col in required_columns
        if col not in df.columns
    ]

    if missing:
        raise ValueError(
            f"{sample_dir}: missing columns {missing}"
        )

    for _, row in df.iterrows():

        # Option B: identify by factor and VF protein
        key = (
            row["Virulence Factor"],
            row["VF Protein"]
        )

        if key not in vir_variants:
            vir_variants[key] = {}

        if sample_dir.name not in vir_variants[key]:
            vir_variants[key][sample_dir.name] = []

        vir_variants[key][sample_dir.name].append({
            "Coverage(%)": row["Coverage(%)"],
            "Depth": row["Depth"]
        })


# Exclude findings shared by all samples
vir_non_shared = {
    key: sample_data
    for key, sample_data in vir_variants.items()
    if len(sample_data) != len(samples)
}


vir_non_shared_calls = sum(
    len(hits)
    for sample_data in vir_non_shared.values()
    for hits in sample_data.values()
)


# Apply virulence thresholds
vir_retained = {}

for key, sample_data in vir_non_shared.items():

    filtered_samples = {}

    for sample_name, hits in sample_data.items():

        passing_hits = []

        for hit in hits:

            # MGAP coverage may be "100.00%"
            coverage = pd.to_numeric(
                str(hit["Coverage(%)"]).replace("%", ""),
                errors="coerce"
            )

            depth = pd.to_numeric(
                hit["Depth"],
                errors="coerce"
            )

            if (
                coverage >= config["virulence"]["min_coverage"]
                and depth >= config["virulence"]["min_depth"]
            ):
                passing_hits.append(hit)

        if passing_hits:
            filtered_samples[sample_name] = passing_hits

    if filtered_samples:
        vir_retained[key] = filtered_samples


vir_final_calls = sum(
    len(hits)
    for sample_data in vir_retained.values()
    for hits in sample_data.values()
)


add_summary(
    "Virulence",
    vir_original,
    vir_non_shared_calls,
    vir_final_calls,
    len(vir_retained),
    (
        f"Coverage >= {config['virulence']['min_coverage']}% "
        f"and Depth >= {config['virulence']['min_depth']}"
    )
)


# ==========================================
# Save summary
# ==========================================

summary = pd.DataFrame(results)

summary = summary.sort_values(
    "Variant_type"
).reset_index(drop=True)

output.parent.mkdir(
    parents=True,
    exist_ok=True
)

summary.to_csv(
    output,
    index=False
)

print("\nVariant comparison and filtering summary:\n")

print(
    summary.to_string(index=False)
)

print(f"\nSaved to {output}")
