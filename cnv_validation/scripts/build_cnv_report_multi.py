import csv
import sys
from collections import defaultdict
from pathlib import Path

manifest = Path(sys.argv[1])
coverage_dir = Path(sys.argv[2])
repeat_file = Path(sys.argv[3])
status_file = Path(sys.argv[4])
output = Path(sys.argv[5])

def read_csv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f))

candidates = read_csv(manifest)
repeats = read_csv(repeat_file)
statuses = read_csv(status_file)

status_by_sample = {
    row["Sample"]: row
    for row in statuses
}

# Read all sample-specific coverage evidence.
samples = sorted({row["Sample"] for row in candidates})

coverage = {}

for sample in samples:

    if sample not in status_by_sample:
        raise ValueError(
            f"No input status recorded for {sample}"
        )

    status = status_by_sample[sample]["Input_status"]

    if status == "SKIPPED":
        print(f"SKIPPED: {sample} — no coverage expected")
        continue

    if status != "READY":
        raise ValueError(
            f"Unexpected status for {sample}: {status}"
        )

    path = coverage_dir / f"{sample}.csv"

    if not path.is_file():
        raise FileNotFoundError(
            f"Coverage missing for READY sample: {sample}"
        )

    for row in read_csv(path):
        key = (row["Candidate_ID"], sample)
        coverage[key] = row

# Group MGAP predictions and repeat alignments by CNV.
by_cnv = defaultdict(list)
by_repeat = defaultdict(list)

for row in candidates:
    by_cnv[row["Candidate_ID"]].append(row)

for row in repeats:
    by_repeat[row["Candidate_ID"]].append(row)


report = [
    "# CNV Validation Report",
    "",
    "## Method",
    "",
    "- Original predictions: MGAP CNV comparison.",
    "- Independent evidence: BWA-MEM alignments and samtools depth.",
    "- Base quality threshold: 20.",
    "- Mapping quality thresholds: 0, 1 and 30.",
    "- Coverage compared with equal-length upstream and downstream regions.",
    "- Repeat analysis: candidate reference sequences aligned back to the complete reference using BWA-MEM.",
    "- Breakpoint support has not yet been systematically evaluated.",
    "",
]

# ==========================================
# Sample-processing summary
# ==========================================

report.extend([
    "## Sample processing summary",
    "",
    "| Sample | Input source | Analysis status |",
    "|---|---|---|",
])

for sample in samples:

    status = status_by_sample[sample]["Input_status"]
    source = status_by_sample[sample]["Input_details"]

    measured = sum(
        (candidate_id, sample) in coverage
        for candidate_id in by_cnv
    )

    if status == "SKIPPED":
        analysis = "Not analyzed — missing input"

    elif measured == len(by_cnv):
        analysis = f"Analyzed — {measured} CNVs"

    else:
        raise ValueError(
            f"Incomplete CNV coverage for {sample}: "
            f"{measured}/{len(by_cnv)}"
        )

    report.append(
        f"| {sample} | {source} | {analysis} |"
    )

report.append("")

for candidate_id, entries in sorted(by_cnv.items()):

    first = entries[0]

    chrom = first["Chromosome"]
    start = first["Start"]
    end = first["End"]

    called = [
        r["Sample"]
        for r in entries
        if r["MGAP_called"] == "Yes"
    ]

    not_called = [
        r["Sample"]
        for r in entries
        if r["MGAP_called"] == "No"
    ]

    report.extend([
        f"## {candidate_id}",
        "",
        f"**Region:** {chrom}:{start}-{end}",
        "",
        f"**MGAP variant type:** {first['Type']}",
        "",
        f"**MGAP-positive samples:** {', '.join(called)}",
        "",
        f"**Not reported by MGAP:** {', '.join(not_called)}",
        "",
        "### Coverage evidence",
        "",
        "| Sample | MGAP called | MGAP normalized RD | Ratio MAPQ 0 | Ratio MAPQ 30 | MAPQ 0 depth contribution |",
        "|---|---|---:|---:|---:|---:|",
    ])

    ratios = []
    mapq0_fractions = []
    flank_warnings = []
    analyzed_called = []
    analyzed_not_called = []

    for entry in entries:
        sample = entry["Sample"]

        # A skipped sample has no independent coverage evidence.
        if (candidate_id, sample) not in coverage:
            report.append(
                f"| {sample} | {entry['MGAP_called']} | "
                f"{entry['MGAP_normalized_RD'] or 'Not reported'} | "
                "Not analyzed | Not analyzed | Not analyzed |"
            )
            continue

        row = coverage[(candidate_id, sample)]

        ratio0 = float(row["Ratio_MAPQ0"])
        ratio30 = float(row["Ratio_MAPQ30"])
        fraction = float(row["MAPQ0_fraction"])

        if entry["MGAP_called"] == "Yes":
            analyzed_called.append(sample)
        else:
            analyzed_not_called.append(sample)

        ratios.append(ratio30)
        mapq0_fractions.append(fraction)

        upstream = float(row["Upstream_MAPQ30"])
        downstream = float(row["Downstream_MAPQ30"])

        lower = min(upstream, downstream)
        higher = max(upstream, downstream)

        if lower == 0 or higher / lower > 2:
            flank_warnings.append(
                f"- {sample}: Uneven flanking coverage "
                f"(upstream={upstream:.2f}, "
                f"downstream={downstream:.2f}). "
                "The local CNV ratio may be unreliable."
            )

        report.append(
            f"| {sample} | {entry['MGAP_called']} | "
            f"{entry['MGAP_normalized_RD'] or 'Not reported'} | "
            f"{ratio0:.3f} | {ratio30:.3f} | "
            f"{fraction * 100:.1f}% |"
        )

    if flank_warnings:
        report.extend([
            "",
            "### Flanking coverage QC warnings",
            "",
        ])
        report.extend(flank_warnings)

    report.extend([
        "",
        "### Reference-repeat evidence",
        "",
    ])

    matches = by_repeat[candidate_id]

    # Distinguish exact full-length matches from other alignments.
    cnv_length = int(end) - int(start) + 1
    exact_matches = 0

    for match in matches:

        cigar = match["CIGAR"]
        nm = match["NM"]

        if cigar == f"{cnv_length}M" and nm == "0":
            match_type = "Full-length exact match"
            exact_matches += 1
        else:
            match_type = "Non-exact or complex alignment"

        report.append(
            f"- {match['Match_chromosome']}:"
            f"{match['Match_position']} "
            f"({match_type}; "
            f"CIGAR={cigar}, "
            f"NM={nm}, "
            f"MAPQ={match['MAPQ']})"
        )

    report.append(
        f"- **Full-length exact reference matches:** "
        f"{exact_matches}"
    )

    # Identify CNV candidates overlapping this region.
    overlapping_candidates = []

    for other_id, other_entries in by_cnv.items():

        if other_id == candidate_id:
            continue

        other = other_entries[0]

        if other["Chromosome"] != chrom:
            continue

        other_start = int(other["Start"])
        other_end = int(other["End"])

        if int(start) <= other_end and other_start <= int(end):
            overlapping_candidates.append(other_id)

    if overlapping_candidates:
        report.extend([
            "",
            "### Overlapping CNV candidates",
            "",
            "**QC warning:** This candidate overlaps "
            + ", ".join(sorted(overlapping_candidates))
            + ". These predictions may represent related "
              "biological events. Local coverage normalization "
              "may also be affected.",
            "",
        ])

    # No samples were analyzed for this candidate.
    if not ratios:
        report.extend([
            "",
            "### Automated interpretation",
            "",
            "**Coverage evidence:** Not available.",
            "",
            "**MGAP agreement:** Cannot be evaluated because "
            "no samples have independent coverage measurements.",
            "",
            "**Biological validation:** Not assessed — missing input data. "
            "The MGAP predictions have not been independently validated.",
            "",
            "---",
            "",
        ])
        continue
    spread = max(ratios) - min(ratios)

    report.extend([
        "",
        "### Automated interpretation",
        "",
        f"- MAPQ 30 coverage-ratio range: {min(ratios):.3f}–{max(ratios):.3f}.",
        f"- Difference between highest and lowest ratio: {spread:.3f}.",
        f"- MAPQ 0 contribution range: "
        f"{min(mapq0_fractions) * 100:.1f}–"
        f"{max(mapq0_fractions) * 100:.1f}%.",
        f"- Reported reference alignments: {len(matches)}.",
        "",
    ])

    if  analyzed_called and analyzed_not_called and spread <= 0.05:
        report.append(
            "**MGAP agreement:** The MAPQ 30 coverage ratios "
            "are similar across all samples, including those "
            "without an MGAP call. The sample-specific MGAP "
            "prediction is not corroborated by this depth comparison."
        )
    else:
        report.append(
            "**MGAP agreement:** Coverage patterns require "
            "further investigation before determining agreement."
        )

    report.extend([
        "",
        "**Biological validation:** Inconclusive. "
        "Reference-repeat alignments and mapping-quality "
        "patterns must be considered. Read-depth evidence "
        "alone does not confirm or exclude this CNV. "
        "Breakpoint and additional structural evidence "
        "remain to be evaluated.",
        "",
        "---",
        "",
    ])

output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("\n".join(report))

print(f"Candidates reported: {len(by_cnv)}")
print(f"Samples compared: {len(samples)}")
print(f"Report saved: {output}")
