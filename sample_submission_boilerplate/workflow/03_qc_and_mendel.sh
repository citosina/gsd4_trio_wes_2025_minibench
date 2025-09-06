#!/usr/bin/env bash
set -euo pipefail

# 03_qc_and_mendel.sh (PLINK-only)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VCF="${ROOT}/results/vcf/trio.vcf.gz"
OUTDIR="${ROOT}/results/qc"
REF="${ROOT}/data/ref/GRCh38_chr3.fa"
mkdir -p "${OUTDIR}"

command -v plink >/dev/null 2>&1 || { echo "ERROR: plink not found in PATH."; exit 1; }
[[ -s "${VCF}" ]] || { echo "ERROR: Missing VCF ${VCF}"; exit 1; }
[[ -s "${REF}" ]] || { echo "ERROR: Missing reference ${REF}"; exit 1; }

# QC summary for SNP/indel counts + Ti/Tv used by 05_evaluate.py
if command -v bcftools >/dev/null 2>&1; then
  echo "[INFO] bcftools stats → variants.bcftools.stats.txt"
  bcftools stats -F "${REF}" -s - "${VCF}" > "${OUTDIR}/variants.bcftools.stats.txt"
else
  echo "[WARN] bcftools not found; Ti/Tv and counts will be unavailable to 05_evaluate.py"
  : > "${OUTDIR}/variants.bcftools.stats.txt"
fi

echo "[INFO] Using PLINK 1.9 for Mendelian check"
# Step A: import VCF (FID=IID=sampleID)
plink \
  --vcf "${VCF}" \
  --double-id \
  --allow-extra-chr \
  --make-bed \
  --out "${OUTDIR}/trio"

# Step B: force all into single family (FAM1)
cat > "${OUTDIR}/ids.update" <<'EOF'
father  father   FAM1  father
mother  mother   FAM1  mother
proband proband  FAM1  proband
EOF

plink \
  --bfile "${OUTDIR}/trio" \
  --update-ids "${OUTDIR}/ids.update" \
  --make-bed \
  --out "${OUTDIR}/trio_fam1"

# Step C: set parents & sex
cat > "${OUTDIR}/parents.update" <<'EOF'
FAM1  father   0       0
FAM1  mother   0       0
FAM1  proband  father  mother
EOF

cat > "${OUTDIR}/sex.update" <<'EOF'
FAM1  father  1
FAM1  mother  2
FAM1  proband 1
EOF

plink \
  --bfile "${OUTDIR}/trio_fam1" \
  --update-parents "${OUTDIR}/parents.update" \
  --update-sex "${OUTDIR}/sex.update" \
  --mendel \
  --out "${OUTDIR}/mendel_plink"

# Normalize: one line per violation to mendel.txt
if [[ -s "${OUTDIR}/mendel_plink.lmendel" ]]; then
  awk 'NR==1{next} {print $1"\t"$2"\t.\t.\t.\t.\t.\t.\tVIOLATION"}' \
    "${OUTDIR}/mendel_plink.lmendel" > "${OUTDIR}/mendel.txt"
else
  echo -e "#CHROM\tPOS\tID\tREF\tALT\tFATHER_GT\tMOTHER_GT\tCHILD_GT\tSTATUS" > "${OUTDIR}/mendel.txt"
fi

echo "[DONE] QC + Mendelian via PLINK → ${OUTDIR}"

