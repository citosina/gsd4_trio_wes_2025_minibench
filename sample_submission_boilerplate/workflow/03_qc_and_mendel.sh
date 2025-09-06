#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VCF="${ROOT}/results/vcf/trio.vcf.gz"
OUTDIR="${ROOT}/results/qc"
mkdir -p "${OUTDIR}"

[[ -s "${VCF}" ]] || { echo "Missing VCF ${VCF}"; exit 1; }

# PED: FamilyID, SampleID, FatherID, MotherID, Sex(1=M,2=F), Phenotype(1=unaff,2=aff,0=unk)
# Asumimos proband afectado, padre/madre no; puedes ajustar si tienes metadatos distintos.
cat > "${OUTDIR}/trio.ped" <<EOF
FAM1  father   0       0       1       1
FAM1  mother   0       0       2       1
FAM1  proband  father  mother  1       2
EOF

# bcftools stats
bcftools stats -F "${ROOT}/data/ref/GRCh38_chr3.fa" -s - "${VCF}" > "${OUTDIR}/variants.bcftools.stats.txt"

# mendelian violations (requiere IDs en VCF iguales a SM del BAM: father,mother,proband)
bcftools +mendelian "${VCF}" -Ped "${OUTDIR}/trio.ped" > "${OUTDIR}/mendel.txt"

echo "QC in ${OUTDIR}"

