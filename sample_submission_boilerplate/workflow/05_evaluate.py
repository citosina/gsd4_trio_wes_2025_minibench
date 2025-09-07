#!/usr/bin/env python3
import json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VCF = ROOT / "results" / "vcf" / "trio.vcf.gz"
STATS = ROOT / "results" / "qc" / "variants.bcftools.stats.txt"
MENDEL = ROOT / "results" / "qc" / "mendel.txt"

def parse_stats(p: Path):
    snps = indels = None
    titv = None
    if not p.exists():
        return snps, indels, titv

    with p.open() as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")

            # SN lines (tab-delimited)
            if parts and parts[0] == "SN" and len(parts) >= 4:
                key = parts[2].strip().lower().rstrip(":")  # <-- strip trailing colon
                if key == "number of snps":
                    try: snps = int(parts[3])
                    except: pass
                elif key == "number of indels":
                    try: indels = int(parts[3])
                    except: pass
                elif key == "ti/tv ratio" and titv is None:
                    try: titv = float(parts[3])
                    except: pass

            # Modern bcftools Ti/Tv from TSTV
            if parts and parts[0] == "TSTV" and len(parts) >= 5 and titv is None:
                try: titv = float(parts[4])
                except: pass

            # Regex fallbacks (robust to odd spacing)
            if snps is None:
                m = re.search(r'\bnumber of SNPs\b[:\s]+(\d+)', line, re.I)
                if m: snps = int(m.group(1))
            if indels is None:
                m = re.search(r'\bnumber of indels\b[:\s]+(\d+)', line, re.I)
                if m: indels = int(m.group(1))
            if titv is None and line.startswith("TSTV"):
                m = re.search(r'^TSTV\s+\S+\s+\S+\s+\S+\s+([0-9.]+)', line)
                if m: titv = float(m.group(1))
    return snps, indels, titv

def count_mendel(p: Path):
    if not p.exists():
        return 0
    n = 0
    with p.open() as fh:
        for line in fh:
            if line.strip() and not line.startswith("#"):
                n += 1
    return n

def main():
    snps, indels, titv = parse_stats(STATS)
    mendel = count_mendel(MENDEL)
    print(json.dumps({
        "vcf": str(VCF),
        "counts": {"snps": snps, "indels": indels},
        "titv": titv,
        "mendel_violations": mendel
    }, indent=2))

if __name__ == "__main__":
    main()

