import csv
import subprocess
import sys
import tempfile
from pathlib import Path


reference = Path(sys.argv[1])
manifest = Path(sys.argv[2])
output = Path(sys.argv[3])

# Read each candidate once, not once per sample.
candidates = {}

with manifest.open(newline="") as f:
    for row in csv.DictReader(f):
        candidate_id = row["Candidate_ID"]

        candidates[candidate_id] = (
            row["Chromosome"],
            int(row["Start"]),
            int(row["End"])
        )

results = []

with tempfile.TemporaryDirectory() as temporary_dir:
    for candidate_id, (chrom, start, end) in candidates.items():

        region = f"{chrom}:{start}-{end}"
        fasta = Path(temporary_dir) / f"{candidate_id}.fa"

        # Extract the candidate sequence.
        extracted = subprocess.run(
            ["samtools", "faidx", str(reference), region],
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )

        fasta.write_text(extracted.stdout)

        # Align it against the complete reference.
        alignment = subprocess.run(
            ["bwa", "mem", "-a", str(reference), str(fasta)],
            stdout=subprocess.PIPE,
            text=True,
            check=True
        )

        for line in alignment.stdout.splitlines():

            if line.startswith("@"):
                continue

            fields = line.split("\t")

            flag = int(fields[1])

            # Skip unmapped alignments.
            if flag & 4:
                continue

            nm = ""

            for field in fields[11:]:
                if field.startswith("NM:i:"):
                    nm = field.split(":")[-1]

            results.append({
                "Candidate_ID": candidate_id,
                "Candidate_region": region,
                "Match_chromosome": fields[2],
                "Match_position": fields[3],
                "MAPQ": fields[4],
                "CIGAR": fields[5],
                "NM": nm,
                "Alignment_flag": flag
            })

        print(f"{candidate_id}: {region} analyzed")

output.parent.mkdir(parents=True, exist_ok=True)

columns = [
    "Candidate_ID",
    "Candidate_region",
    "Match_chromosome",
    "Match_position",
    "MAPQ",
    "CIGAR",
    "NM",
    "Alignment_flag"
]

with output.open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=columns)
    writer.writeheader()
    writer.writerows(results)

print(f"Total reference alignments: {len(results)}")
print(f"Evidence saved: {output}")
