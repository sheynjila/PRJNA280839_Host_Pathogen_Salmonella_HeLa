#!/bin/bash
#SBATCH --job-name=sh06_quant
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00

###############################################################################
# 6_Split_And_Quantify_Host_Pathogen.sh
#
# NARRATIVE (the second half of this pipeline's core host/pathogen design):
# Script 5 produced ONE BAM file containing a mixture of host- and
# pathogen-mapped reads. This script separates that mixture into a host
# expression profile and a pathogen expression profile WITHOUT physically
# splitting the BAM: featureCounts is run TWICE on the SAME combined BAM,
# once against the host-only GTF, once against the pathogen-only GTF.
# featureCounts only counts a read if it overlaps a feature in the GTF
# it's given, so a host-only GTF simply never assigns a pathogen-mapped
# read to a gene (and vice versa).
#
# UNLIKE the companion PRJNA1091619 pipeline's single-contig Mtb pathogen
# genome, Salmonella Typhimurium SL1344's reference has FOUR pathogen
# contigs: the main chromosome (NC_016810.1) plus three plasmids
# (NC_017718.1, NC_017719.1, NC_017720.1) -- all four, not just one, must
# be excluded from "host" and counted as "pathogen" in the read-fraction
# diagnostic below.
###############################################################################
set -e
set -o pipefail

module purge
module load SAMtools/1.19.2-GCC-13.2.0

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_stepbystep"
REFDIR="/scratch/$(whoami)/PRJNA280839_salmonella_ref"
mkdir -p "${WORKDIR}"/{counts,qc_reports}
cd "${WORKDIR}"

SRR="SRR2936892"
BAM="alignments/${SRR}_Aligned.sortedByCoord.out.bam"

# ------------------------------------------------------------------------------
# Diagnostic: host vs pathogen mapped-read fraction. Pathogen contigs are
# the SL1344 chromosome + its three plasmids -- listed explicitly, not
# assumed to be "everything that isn't a human chromosome" the way the
# single-contig Mtb pipeline could get away with.
# ------------------------------------------------------------------------------
samtools idxstats "${BAM}" > "qc_reports/${SRR}_idxstats.txt"

PATHOGEN_CONTIGS="NC_016810.1|NC_017718.1|NC_017719.1|NC_017720.1"
PATHOGEN_MAPPED=$(awk -v pat="${PATHOGEN_CONTIGS}" '$1 ~ "^(" pat ")$" {sum += $3} END {print sum+0}' "qc_reports/${SRR}_idxstats.txt")
HOST_MAPPED=$(awk -v pat="${PATHOGEN_CONTIGS}" '$1 !~ "^(" pat ")$" && $1 != "*" {sum += $3} END {print sum+0}' "qc_reports/${SRR}_idxstats.txt")
TOTAL_MAPPED=$((PATHOGEN_MAPPED + HOST_MAPPED))

{
    echo "# Host vs pathogen mapped-read fraction: ${SRR}"
    echo "Host-mapped reads:     ${HOST_MAPPED}"
    echo "Pathogen-mapped reads: ${PATHOGEN_MAPPED}"
    echo "Total mapped reads:    ${TOTAL_MAPPED}"
    if [ "${TOTAL_MAPPED}" -gt 0 ]; then
        awk -v p="${PATHOGEN_MAPPED}" -v t="${TOTAL_MAPPED}" 'BEGIN { printf "Pathogen fraction:      %.4f%%\n", (p/t)*100 }'
    fi
} | tee "qc_reports/${SRR}_host_pathogen_fraction.txt"

# ------------------------------------------------------------------------------
# Quantification, HOST side. -s 0 (unstranded): confirm against your own
# data's actual strandedness before trusting this on a different dataset.
# ------------------------------------------------------------------------------
module purge
module load Subread/2.0.6-GCC-12.3.0

featureCounts -T 4 -s 0 \
    -a "${REFDIR}/host_GRCh38.gtf" \
    -o "counts/${SRR}_host_counts.txt" \
    "${BAM}"

# ------------------------------------------------------------------------------
# Quantification, PATHOGEN side: same BAM, pathogen-only GTF. This is
# where the paper's headline finding lives -- the small RNA "PinT" is a
# Salmonella sRNA feature in this GTF, not a CDS.
# ------------------------------------------------------------------------------
featureCounts -T 4 -s 0 \
    -a "${REFDIR}/pathogen_SL1344.gtf" \
    -o "counts/${SRR}_pathogen_counts.txt" \
    "${BAM}"

echo "OK: host and pathogen count matrices written to counts/."
echo "    counts/${SRR}_host_counts.txt.summary and _pathogen_counts.txt.summary"
echo "    report each side's featureCounts assignment rate."
