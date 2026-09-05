#!/bin/bash
#SBATCH --job-name=sh02_fastqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=00:45:00

###############################################################################
# 2_Interpret_FastQC_Diagnostics.sh
#
# NARRATIVE:
# Runs FastQC on the single-end raw FASTQ from Script 1, then parses its
# summary.txt into PASS/WARN/FAIL counts -- the evidence Script 4 uses to
# decide whether to enable adapter/quality trimming, not a blind default.
###############################################################################
set -e
set -o pipefail

module purge
module load FastQC/0.11.9-Java-11

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep"
mkdir -p "${WORKDIR}/qc_before"
cd "${WORKDIR}"

SRR="SRR2936892"

echo "Running FastQC on ${SRR} (single-end)..."
fastqc "raw_reads/${SRR}.fastq" -o qc_before -t 4 --extract

SUMMARY="qc_before/${SRR}_fastqc/summary.txt"
echo "---- ${SRR}: FastQC verdicts ----"
cat "${SUMMARY}"

FAIL_COUNT=$(grep -c "^FAIL" "${SUMMARY}" || true)
WARN_COUNT=$(grep -c "^WARN" "${SUMMARY}" || true)
echo "Synthesis: ${FAIL_COUNT} FAIL, ${WARN_COUNT} WARN."
