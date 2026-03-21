#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/mutmap_out_$(date +%Y%m%d_%H%M%S)"

mutmap -o "${OUT_DIR}" \
       -n 20 \
       -w 100 \
       -s 20 \
       -r "${SCRIPT_DIR}/mutmap_ref.fasta" \
       -c "${SCRIPT_DIR}/mutmap_cultivar.1.fastq.gz,${SCRIPT_DIR}/mutmap_cultivar.2.fastq.gz" \
       -b "${SCRIPT_DIR}/mutmap_bulk.1.fastq.gz,${SCRIPT_DIR}/mutmap_bulk.2.fastq.gz"

echo "MutMap test completed. Output: ${OUT_DIR}"
