#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ="${ROOT}/data/fastq"
REF="${ROOT}/data/ref/GRCh38_chr3.fa"
OUT="${ROOT}/results/bam"
TMP="${ROOT}/results/tmp"
mkdir -p "${OUT}" "${TMP}"

# Map SRR -> sample label (RGSM) without associative arrays (Bash 3 compatible)
sm_label() {
  case "$1" in
    SRR27997290) echo "father" ;;
    SRR27997291) echo "mother" ;;
    SRR27997292) echo "proband" ;;
    *) echo "sample_$1" ;;
  esac
}

for SID in SRR27997290 SRR27997291 SRR27997292; do
  R1="${FASTQ}/${SID}_1.mini150k_s100.fastq.gz"
  R2="${FASTQ}/${SID}_2.mini150k_s100.fastq.gz"
  if [[ ! -s "$R1" || ! -s "$R2" ]]; then
    echo "ERROR: Missing FASTQs for ${SID} at:"
    echo "  $R1"
    echo "  $R2"
    exit 1
  fi

  RGID="${SID}"
  RGPU="${SID}.PU1"
  RGLB="WES_chr3"
  RGSM="$(sm_label "${SID}")"
  RGPL="ILLUMINA"

  BAM="${OUT}/${SID}.sorted.bam"
  MARKED="${OUT}/${SID}.mkdup.bam"

  echo "[${SID}] bwa-mem2 mem → sort"
  bwa-mem2 mem -t 8 \
    -R "@RG\tID:${RGID}\tSM:${RGSM}\tPL:${RGPL}\tLB:${RGLB}\tPU:${RGPU}" \
    "${REF}" "${R1}" "${R2}" \
    | samtools sort -@4 -o "${BAM}" -

  echo "[${SID}] index + markdup"
  samtools index -@2 "${BAM}"
  picard MarkDuplicates \
    I="${BAM}" O="${MARKED}" \
    M="${OUT}/${SID}.mkdup.metrics.txt" \
    CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT

  rm -f "${BAM}" "${BAM}.bai"
  echo "[${SID}] done → ${MARKED}"
done

echo "All BAMs ready in ${OUT}"

