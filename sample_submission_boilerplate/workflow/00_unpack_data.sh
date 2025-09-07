#!/usr/bin/env bash
set -euo pipefail

# 00_unpack_data.sh
# Usage:
#   bash 00_unpack_data.sh /path/to/GSD4_trio_DataFiles.zip
#   # or, if already extracted:
#   bash 00_unpack_data.sh /path/to/GSD4_trio   # contains fastq/ and ref/

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${HERE%/workflow}"
DATA_DIR="${ROOT}/data"
DEST_FASTQ="${DATA_DIR}/fastq"
DEST_REF="${DATA_DIR}/ref"

SRC="${1:-}"
if [[ -z "${SRC}" ]]; then
  echo "ERROR: please pass the path to GSD4_trio_DataFiles.zip (or an extracted folder with fastq/ and ref/)" >&2
  exit 1
fi

mkdir -p "${DEST_FASTQ}" "${DEST_REF}"

WORKDIR=""
CLEANUP=""
ROOT_IN_ZIP=""

if [[ -f "${SRC}" && "${SRC}" == *.zip ]]; then
  # Extract to a temp dir
  WORKDIR="$(mktemp -d)"
  CLEANUP="yes"
  echo "[INFO] Extracting '${SRC}' to: ${WORKDIR}"
  unzip -q "${SRC}" -d "${WORKDIR}"

  # Try common layout: files at top-level
  if [[ -d "${WORKDIR}/fastq" && -d "${WORKDIR}/ref" ]]; then
    ROOT_IN_ZIP="${WORKDIR}"
  else
    # Find a directory that has both fastq/ and ref/
    ROOT_IN_ZIP="$(find "${WORKDIR}" -maxdepth 2 -type d -name fastq -exec dirname {} \; | head -n1 || true)"
    if [[ -z "${ROOT_IN_ZIP}" || ! -d "${ROOT_IN_ZIP}/fastq" || ! -d "${ROOT_IN_ZIP}/ref" ]]; then
      echo "ERROR: Could not locate fastq/ and ref/ inside the zip." >&2
      exit 1
    fi
  fi
  SRC_ROOT="${ROOT_IN_ZIP}"
else
  # Assume folder containing fastq/ and ref/
  if [[ -d "${SRC}/fastq" && -d "${SRC}/ref" ]]; then
    SRC_ROOT="${SRC}"
  else
    echo "ERROR: '${SRC}' is not a zip and does not contain fastq/ and ref/." >&2
    exit 1
  fi
fi

echo "[INFO] Copying FASTQs and reference to ${DATA_DIR}/"
cp -a "${SRC_ROOT}/fastq/." "${DEST_FASTQ}/"
cp -a "${SRC_ROOT}/ref/."   "${DEST_REF}/"

# Verify checksums if present at source; otherwise compute new ones in destination
if [[ -f "${SRC_ROOT}/CHECKSUMS.sha256" ]]; then
  echo "[INFO] Verifying checksums from source..."
  # Recompute in destination so path prefixes match
  (cd "${DATA_DIR}" && find fastq ref -type f -print0 | xargs -0 sha256sum > CHECKSUMS.sha256)
  (cd "${DATA_DIR}" && shasum -a 256 -c CHECKSUMS.sha256)
  echo "[INFO] Checksums OK."
else
  echo "[WARN] No CHECKSUMS.sha256 found in source; generating fresh checksums in ${DATA_DIR}"
  (cd "${DATA_DIR}" && find fastq ref -type f -print0 | xargs -0 sha256sum > CHECKSUMS.sha256)
fi

echo "[INFO] Done. Data is ready in: ${DATA_DIR}"
# cleanup temp dir
if [[ -n "${CLEANUP}" ]]; then
  rm -rf "${WORKDIR}"
fi
