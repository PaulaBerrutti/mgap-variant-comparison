import csv
import sys
from pathlib import Path


input_csv = Path(sys.argv[1])
output_csv = Path(sys.argv[2])

samples = sys.argv[3:]

if not samples:
    raise ValueError("Provide at least one sample name.")

output_csv.parent.mkdir(parents=True, exist_ok=True)

with input_csv.open(newline="") as infile:
    candidates = list(csv.DictReader(infile))

with output_csv.open("w", newline="") as outfile:
    writer = csv.writer(outfile)

    writer.writerow([
        "Candidate_ID",
        "Chromosome",
        "Start",
        "End",
        "Type",
        "Sample",
        "MGAP_called",
        "MGAP_normalized_RD",
    ])

    for number, candidate in enumerate(candidates, start=1):

        candidate_id = f"CNV{number:03d}"

        called_samples = {
            sample.strip()
            for sample in candidate["Present_in_samples"].split(",")
        }

        rd_by_sample = {}

        for entry in candidate["Normalized_RD_by_sample"].split(";"):
            if entry.strip():
                sample, rd = entry.split(":", 1)
                rd_by_sample[sample.strip()] = rd.strip()

        if not called_samples.issubset(set(samples)):
            raise ValueError(
                f"Unknown sample in {candidate_id}: {called_samples - set(samples)}"
            )

        for sample in samples:
            writer.writerow([
                candidate_id,
                candidate["Chromosome"],
                candidate["Start"],
                candidate["End"],
                candidate["Type"],
                sample,
                "Yes" if sample in called_samples else "No",
                rd_by_sample.get(sample, ""),
            ])

print(f"Imported {len(candidates)} CNV candidates.")
print(f"Created {len(candidates) * len(samples)} candidate/sample rows.")
print(f"Output: {output_csv}")
