import csv
import subprocess
import sys
from pathlib import Path

sample = sys.argv[1]
bam = Path(sys.argv[2])
manifest = Path(sys.argv[3])
output = Path(sys.argv[4])

with manifest.open(newline="") as f:
    candidates = [
        row for row in csv.DictReader(f)
        if row["Sample"] == sample
    ]

if not candidates:
    raise ValueError(f"No CNV candidates found for {sample}")


def measure_depth(chrom, start, end, mapq):
    """Calculate mean depth in three equally sized regions."""

    length = end - start + 1
    region_start = start - length
    region_end = end + length

    if region_start < 1:
        raise ValueError("Upstream region extends before position 1")

    command = [
        "samtools", "depth",
        "-a", "-s",
        "-q", "20",
        "-Q", str(mapq),
        "-r", f"{chrom}:{region_start}-{region_end}",
        str(bam)
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=True
    )

    upstream = 0
    candidate = 0
    downstream = 0
    count = 0

    for line in result.stdout.splitlines():
        ref, pos, depth = line.split()
        pos = int(pos)
        depth = int(depth)
        count += 1

        if pos < start:
            upstream += depth
        elif pos <= end:
            candidate += depth
        else:
            downstream += depth

    if count != 3 * length:
        raise ValueError(
            f"Expected {3 * length} positions, found {count}"
        )

    up_mean = upstream / length
    cnv_mean = candidate / length
    down_mean = downstream / length
    flank_mean = (upstream + downstream) / (2 * length)

    return {
        "upstream": up_mean,
        "candidate": cnv_mean,
        "downstream": down_mean,
        "ratio": (
            cnv_mean / flank_mean
            if flank_mean > 0 else None
        )
    }


output.parent.mkdir(parents=True, exist_ok=True)

rows = []

for candidate in candidates:
    chrom = candidate["Chromosome"]
    start = int(candidate["Start"])
    end = int(candidate["End"])

    measurements = {}

    for mq in [0, 1, 30]:
        measurements[mq] = measure_depth(
            chrom, start, end, mq
        )

    depth0 = measurements[0]["candidate"]
    depth1 = measurements[1]["candidate"]

    mapq0_fraction = (
        (depth0 - depth1) / depth0
        if depth0 > 0 else None
    )

    row = {
        **candidate,
        "Depth_MAPQ0": depth0,
        "Depth_MAPQ1": depth1,
        "Depth_MAPQ30": measurements[30]["candidate"],
        "Upstream_MAPQ30": measurements[30]["upstream"],
        "Downstream_MAPQ30": measurements[30]["downstream"],
        "Ratio_MAPQ0": measurements[0]["ratio"],
        "Ratio_MAPQ30": measurements[30]["ratio"],
        "MAPQ0_fraction": mapq0_fraction
    }

    rows.append(row)

with output.open("w", newline="") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(rows[0].keys())
    )
    writer.writeheader()
    writer.writerows(rows)

print(f"Sample: {sample}")
print(f"CNVs analyzed: {len(rows)}")

for row in rows:
    print(
        row["Candidate_ID"],
        "MAPQ30 ratio:",
        row["Ratio_MAPQ30"],
        "MAPQ0 fraction:",
        row["MAPQ0_fraction"]
    )

print(f"Saved: {output}")
