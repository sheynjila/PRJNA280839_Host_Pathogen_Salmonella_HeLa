#!/bin/bash
#SBATCH --job-name=sh04_fastp
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

###############################################################################
# 4_Evidence_Based_Fastp_Trimming.sh
#
# NARRATIVE:
# Evidence-based trimming: only enables adapter/quality flags Script 2's
# actual FastQC verdict justifies. Single-end fastp (-i/-o), matching this
# study's real library layout (verified via SRA run metadata).
###############################################################################
set -e
set -o pipefail

module purge
module load FastQC/0.11.9-Java-11
module load fastp/0.23.4-GCC-13.2.0

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep"
mkdir -p "${WORKDIR}"/{trimmed_reads,fastp_reports,qc_after}
cd "${WORKDIR}"

SRR="SRR2936892"
SUMMARY="qc_before/${SRR}_fastqc/summary.txt"

FASTP_FLAGS=("--length_required" "36")
RATIONALE_LINES=("- Minimum length floor of 36bp applied unconditionally.")

if grep -qE "^(WARN|FAIL)\s+Adapter Content" "${SUMMARY}"; then
    FASTP_FLAGS+=("--detect_adapter_for_pe")
    RATIONALE_LINES+=("- Adapter trimming ENABLED: FastQC flagged Adapter Content.")
fi

if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "${SUMMARY}"; then
    FASTP_FLAGS+=("--cut_tail" "--cut_tail_mean_quality" "20")
    RATIONALE_LINES+=("- 3' quality trimming ENABLED: FastQC flagged quality decay.")
fi

echo "Running single-end fastp with flags: ${FASTP_FLAGS[*]}"
fastp -i "raw_reads/${SRR}.fastq" \
      -o "trimmed_reads/${SRR}_clean.fastq" \
      "${FASTP_FLAGS[@]}" --thread 4 \
      --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html"

{
    echo "# Trimming rationale for ${SRR}"
    printf '%s\n' "${RATIONALE_LINES[@]}"
} > "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"

fastqc "trimmed_reads/${SRR}_clean.fastq" -o qc_after -t 4 --extract
