#!/bin/bash
#SBATCH --job-name=sh03_multiqc
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:30:00

###############################################################################
# 3_Batch_MultiQC_Summary.sh
#
# NARRATIVE:
# Aggregates qc_before/ into one dashboard. Re-run after Script 7 (the
# automated all-8-samples run) to see every sample's QC side by side.
###############################################################################
set -e
set -o pipefail

module purge
module load MultiQC/1.28-foss-2024a

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep"
cd "${WORKDIR}"

echo "Aggregating qc_before/ reports with MultiQC..."
multiqc qc_before/ -n multiqc_raw_batch.html -o . --force

echo "OK: Combined report generated as multiqc_raw_batch.html"
