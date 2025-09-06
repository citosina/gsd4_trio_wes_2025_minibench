#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VCF_IN="${ROOT}/results/vcf/trio.vcf.gz"
ANN_DIR="${ROOT}/results/ann"
mkdir -p "${ANN_DIR}"

[[ -s "${VCF_IN}" ]] || { echo "Missing VCF ${VCF_IN}"; exit 1; }

SNPEFF_DB="GRCh38.86"

snpEff -canon -noStats "${SNPEFF_DB}" "${VCF_IN}" | bgzip -c > "${ANN_DIR}/trio.ann.vcf.gz"
bcftools index -f "${ANN_DIR}/trio.ann.vcf.gz"

bcftools view -i 'ANN~"\\|GBE1\\|"' "${ANN_DIR}/trio.ann.vcf.gz" -Oz -o "${ANN_DIR}/GBE1_candidates.vcf.gz"
bcftools index -f "${ANN_DIR}/GBE1_candidates.vcf.gz"

echo "Annotation written to ${ANN_DIR}"

