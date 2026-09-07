#!/bin/bash
#SBATCH --job-name=sh05_align
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=03:00:00

###############################################################################
# 5_Build_Dual_Reference_And_Align.sh
#
# NARRATIVE (the core host-vs-pathogen design decision of this pipeline):
# This is a dual RNA-seq experiment: every FASTQ file is a mixture of reads
# from the Salmonella-infected HeLa cell (host) AND from the Salmonella
# Typhimurium SL1344 bacteria infecting it (pathogen). There is no way to
# know, before alignment, which read came from which organism.
#
# The standard, simplest way to resolve this -- the SAME approach this
# curriculum's companion PRJNA1091619 pipeline uses, which itself named
# THIS STUDY (Westermann et al. 2016) as its methodological precedent -- is
# a CONCATENATED REFERENCE: build ONE STAR index over both genomes glued
# together, and align every read against it in a single pass. STAR then
# assigns each read to whichever genome it actually matches. The actual
# host-vs-pathogen SEPARATION happens downstream, in Script 6, at the
# gene-counting stage, not here.
#
# References used (verified real accessions, not placeholders):
#   Host:     Ensembl GRCh38, release 110.
#   Pathogen: NCBI RefSeq GCF_000210855.2 (assembly ASM21085v2), the
#             Salmonella enterica subsp. enterica serovar Typhimurium
#             str. SL1344 reference genome (chromosome NC_016810.1 plus
#             three plasmids, NC_017718/017719/017720).
#
# MEMORY -- APPLYING A REAL LESSON FROM THE START, NOT WAITING TO REDISCOVER
# IT: this curriculum's companion PRJNA1091619 pipeline was real-hardware
# OOM-killed TWICE building a combined host+pathogen STAR index at
# --mem=36G, and separately discovered Slurm can grant less memory than
# requested even when the request itself is correct (see that pipeline's
# TROUBLESHOOTING.md, Failures 1 and 3). This script starts at --mem=64G
# and includes the same startup diagnostic from the outset, rather than
# repeating that entire discovery process on a new dataset.
#
# --sjdbOverhang 74 = (this study's own real 75bp single-end read length,
# verified via SRA run metadata) - 1.
###############################################################################
set -e
set -o pipefail

module purge
module load STAR/2.7.10b-GCC-11.3.0

echo "[$(date +'%b %d %H:%M:%S')] Slurm-granted memory for this job: ${SLURM_MEM_PER_NODE:-unknown} MB (this script requested 65536 MB via --mem=64G)"
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE}" -lt 65536 ]; then
    echo "[$(date +'%b %d %H:%M:%S')] WARNING: Slurm granted less memory than this script requested." >&2
    echo "  Requested via #SBATCH --mem=64G (65536 MB); actually granted ${SLURM_MEM_PER_NODE} MB." >&2
    echo "  This is a submission/cluster-side issue (a --mem override on the sbatch command line," >&2
    echo "  or a partition/QOS/association memory cap), not a problem with this script's request." >&2
    echo "  Check: sacct -j \${SLURM_JOB_ID} -o JobID,ReqMem,AllocTRES%40 and scontrol show partition <name>." >&2
fi

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep_${RUN_TAG:-orchestrated}"
REFDIR="/scratch/$(whoami)/PRJNA280839_salmonella_ref"
mkdir -p "${WORKDIR}/alignments" "${REFDIR}"
cd "${WORKDIR}"

HOST_FASTA_URL="http://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
HOST_GTF_URL="http://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"
PATHOGEN_FASTA_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.fna.gz"
PATHOGEN_GTF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/210/855/GCF_000210855.2_ASM21085v2/GCF_000210855.2_ASM21085v2_genomic.gtf.gz"

# ------------------------------------------------------------------------------
# Step A: download each reference SEPARATELY first (keeps host and pathogen
# annotation files available on their own for Script 6's two featureCounts
# calls -- see that script for why each species needs its own GTF).
# ------------------------------------------------------------------------------
if [ ! -f "${REFDIR}/host_GRCh38.fa" ]; then
    wget -qO- "${HOST_FASTA_URL}" | gunzip > "${REFDIR}/host_GRCh38.fa"
fi
if [ ! -f "${REFDIR}/host_GRCh38.gtf" ]; then
    wget -qO- "${HOST_GTF_URL}" | gunzip > "${REFDIR}/host_GRCh38.gtf"
fi
if [ ! -f "${REFDIR}/pathogen_SL1344.fa" ]; then
    wget -qO- "${PATHOGEN_FASTA_URL}" | gunzip > "${REFDIR}/pathogen_SL1344.fa"
fi
if [ ! -f "${REFDIR}/pathogen_SL1344.gtf" ]; then
    wget -qO- "${PATHOGEN_GTF_URL}" | gunzip > "${REFDIR}/pathogen_SL1344.gtf"
fi

# ------------------------------------------------------------------------------
# Step B: concatenate into ONE combined reference for STAR's single-pass
# alignment (Step C).
# ------------------------------------------------------------------------------
if [ ! -f "${REFDIR}/combined_host_pathogen.fa" ]; then
    cat "${REFDIR}/host_GRCh38.fa" "${REFDIR}/pathogen_SL1344.fa" > "${REFDIR}/combined_host_pathogen.fa"
fi
if [ ! -f "${REFDIR}/combined_host_pathogen.gtf" ]; then
    cat "${REFDIR}/host_GRCh38.gtf" "${REFDIR}/pathogen_SL1344.gtf" > "${REFDIR}/combined_host_pathogen.gtf"
fi

# ------------------------------------------------------------------------------
# Step C: build the combined STAR index (idempotent -- checks for
# genomeParameters.txt specifically, the file STAR only writes at the very
# end of a SUCCESSFUL build, not just directory non-emptiness. See this
# curriculum's PRJNA1091619 pipeline TROUBLESHOOTING.md, Failure 2, for the
# real stale-partial-index failure this proactively guards against.
# ------------------------------------------------------------------------------
if [ ! -f "${REFDIR}/star_index/genomeParameters.txt" ]; then
    rm -rf "${REFDIR}/star_index"
    mkdir -p "${REFDIR}/star_index"
    STAR --runThreadN 8 --runMode genomeGenerate --genomeDir "${REFDIR}/star_index" \
         --genomeFastaFiles "${REFDIR}/combined_host_pathogen.fa" \
         --sjdbGTFfile "${REFDIR}/combined_host_pathogen.gtf" \
         --sjdbOverhang 74
fi

# ------------------------------------------------------------------------------
# Step D: align the representative sample. Single-end alignment: one file
# passed to --readFilesIn, never two.
# ------------------------------------------------------------------------------
SRR="SRR2936892"

echo "Aligning ${SRR} (single-end) against the combined host+pathogen index..."
STAR --runThreadN 8 --genomeDir "${REFDIR}/star_index" \
     --readFilesIn "trimmed_reads/${SRR}_clean.fastq" \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix "alignments/${SRR}_"
