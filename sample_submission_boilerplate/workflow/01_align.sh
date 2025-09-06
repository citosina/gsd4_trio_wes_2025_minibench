#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FASTQ="${ROOT}/data/fastq"
REF="${ROOT}/data/ref/GRCh38_chr3.fa"
OUT="${ROOT}/results/bam"
TMP="${ROOT}/results/tmp"
mkdir -p "${OUT}" "${TMP}"

# Mapeo SRR -> etiqueta de muestra para SM (padre/madre/hijo)
declare -A SM
SM[SRR27997290]="father"
SM[SRR27997291]="mother"
SM[SRR27997292]="proband"

for SID in SRR27997290 SRR27997291 SRR27997292; do
  R1="${FASTQ}/${SID}_1.mini150k_s100.fastq.gz"
  R2="${FASTQ}/${SID}_2.mini150k_s100.fastq.gz"
  [[ -s "$R1" && -s "$R2" ]] || { echo "Missing FASTQs for ${SID}"; exit 1; }

  RGID="${SID}"
  RGPU="${SID}.PU1"
  RGLB="WES_chr3"
  RGSM="${SM[$SID]}"
  RGPL="ILLUMINA"

  BAM="${OUT}/${SID}.sorted.bam"
  MARKED="${OUT}/${SID}.mkdup.bam"

  echo "[${SID}] bwa-mem2..."
  bwa-mem2 mem -t 8 -R "@RG\tID:${RGID}\tSM:${RGSM}\tPL:${RGPL}\tLB:${RGLB}\tPU:${RGPU}" \
    "${REF}" "${R1}" "${R2}" \
    | samtools sort -@4 -o "${BAM}" -

  echo "[${SID}] markdup..."
  samtools index -@2 "${BAM}"
  picard MarkDuplicates I="${BAM}" O="${MARKED}" M="${OUT}/${SID}.mkdup.metrics.txt" CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  rm -f "${BAM}" "${BAM}.bai"

  echo "[${SID}] done -> ${MARKED}"
done

echo "All BAMs in ${OUT}"

