#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF="${ROOT}/data/ref/GRCh38_chr3.fa"
BAMDIR="${ROOT}/results/bam"
VCFDIR="${ROOT}/results/vcf"
mkdir -p "${VCFDIR}"

BAMS=(
  "${BAMDIR}/SRR27997290.mkdup.bam"
  "${BAMDIR}/SRR27997291.mkdup.bam"
  "${BAMDIR}/SRR27997292.mkdup.bam"
)
for b in "${BAMS[@]}"; do [[ -s "$b" ]] || { echo "Missing BAM: $b"; exit 1; }; done

VCF="${VCFDIR}/trio.vcf.gz"

bcftools mpileup -Ou -f "${REF}" -a DP,AD,ADF,ADR -q 20 -Q 20 "${BAMS[@]}" \
  | bcftools call -mv -Oz -o "${VCF}"

bcftools index -f "${VCF}"
echo "Wrote ${VCF}"

