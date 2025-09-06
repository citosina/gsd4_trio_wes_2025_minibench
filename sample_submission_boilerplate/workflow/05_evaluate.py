#!/usr/bin/env python3
import json, re, sys, pathlib

root = pathlib.Path(__file__).resolve().parents[1]
qc_dir = root / "results" / "qc"
stats_file = qc_dir / "variants.bcftools.stats.txt"
mendel_file = qc_dir / "mendel.txt"
vcf_file = root / "results" / "vcf" / "trio.vcf.gz"

res = {
    "vcf": str(vcf_file),
    "counts": {"snps": 0, "indels": 0},
    "titv": None,
    "mendel_violations": 0
}

# Parse bcftools stats
if stats_file.exists():
    with open(stats_file) as fh:
        for line in fh:
            if line.startswith("SN"):
                if "number of SNPs:" in line:
                    res["counts"]["snps"] = int(line.strip().split("\t")[-1])
                elif "number of indels:" in line:
                    res["counts"]["indels"] = int(line.strip().split("\t")[-1])
                elif "Ti/Tv ratio" in line:
                    try:
                        res["titv"] = float(line.strip().split("\t")[-1])
                    except ValueError:
                        pass

# Parse mendelian violations
if mendel_file.exists():
    with open(mendel_file) as fh:
        res["mendel_violations"] = sum(1 for l in fh if "VIOLATION" in l)

print(json.dumps(res, indent=2))

