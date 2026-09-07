#!/bin/bash
#SBATCH --job-name=sh01_fetch
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00

###############################################################################
# 1_Fetch_And_Verify_Raw_Reads.sh
#
# NARRATIVE:
# This is a dual RNA-seq experiment: HeLa-S3 cells infected with Salmonella
# enterica serovar Typhimurium SL1344, sequenced single-end -- every FASTQ
# is a mixture of human host reads and Salmonella pathogen reads (see
# Script 5's narrative for the combined-reference design this implies).
#
# Representative sample: SRR2936892 (GSM1924778, Wild-type Salmonella
# infection, 4h post-infection, replicate 1). Real BioProject PRJNA280839,
# SRA study SRP057065, verified against SRA run metadata directly (75bp
# single-end, Illumina NextSeq 500).
###############################################################################
set -e
set -o pipefail

module purge
module load SRA-Toolkit/3.0.3-gompi-2022a

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep_${RUN_TAG:-orchestrated}"
mkdir -p "${WORKDIR}"/{raw_reads,checksums,sra_cache}
cd "${WORKDIR}"

SRR="SRR2936892"

echo "Fetching ${SRR} (single-end, WT Salmonella infection 4h, replicate 1)..."
prefetch "${SRR}" --output-directory sra_cache
vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | grep -q "is consistent"
fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --outdir raw_reads --threads 4

md5sum "raw_reads/${SRR}.fastq" > "checksums/${SRR}_raw.md5"

echo "OK: raw_reads/${SRR}.fastq fetched and checksummed."
