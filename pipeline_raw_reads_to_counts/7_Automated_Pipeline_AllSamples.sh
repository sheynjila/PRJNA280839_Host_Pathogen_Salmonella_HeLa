#!/bin/bash
#SBATCH --job-name=sh07_pipeline
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

###############################################################################
# 7_Automated_Pipeline_AllSamples.sh
#
# NARRATIVE:
# Production entry point: loops fetch -> FastQC -> evidence-based fastp ->
# dual-reference alignment -> host/pathogen split-quantification over all 8
# samples in the representative subset (data/metadata_8samples.csv: WT
# Salmonella infection vs mock, at 4h and 24h, 2 replicates each -- all
# single-end, verified via SRA run metadata, so no layout branching is
# needed here, unlike this curriculum's mixed-layout pasilla pipeline).
#
# HARD PREREQUISITE: this script assumes the shared dual-reference STAR
# index at ${REFDIR}/star_index already exists. Run Script 5 once first.
#
# FAULT TOLERANCE, DONE CORRECTLY: `if ! ( subshell ); then ...`, not a
# bare subshell statement -- see PRJNA229998_Airway_Pipeline's README for
# the real `set -e` interaction bug this specific construction avoids.
###############################################################################
set -e
set -o pipefail

echo "[$(date +'%b %d %H:%M:%S')] Slurm-granted memory for this job: ${SLURM_MEM_PER_NODE:-unknown} MB (this script requested 65536 MB via --mem=64G)"
if [ -n "${SLURM_MEM_PER_NODE:-}" ] && [ "${SLURM_MEM_PER_NODE}" -lt 65536 ]; then
    echo "[$(date +'%b %d %H:%M:%S')] WARNING: Slurm granted less memory than this script requested." >&2
    echo "  Requested via #SBATCH --mem=64G (65536 MB); actually granted ${SLURM_MEM_PER_NODE} MB." >&2
    echo "  This is a submission/cluster-side issue (a --mem override on the sbatch command line," >&2
    echo "  or a partition/QOS/association memory cap), not a problem with this script's request." >&2
    echo "  Check: sacct -j \${SLURM_JOB_ID} -o JobID,ReqMem,AllocTRES%40 and scontrol show partition <name>." >&2
fi

WORKDIR="/scratch/$(whoami)/PRJNA280839_salmonella_automated"
REFDIR="/scratch/$(whoami)/PRJNA280839_salmonella_ref"
mkdir -p "${WORKDIR}"/{raw_reads,checksums,qc_before,trimmed_reads,fastp_reports,qc_after,alignments,counts,qc_reports,logs,sra_cache}
cd "${WORKDIR}"

if [ ! -f "${REFDIR}/star_index/genomeParameters.txt" ]; then
    echo "ERROR: ${REFDIR}/star_index is missing or incomplete. Run Script 5 first (and if a star_index/ directory already exists there, delete it first)." >&2
    exit 1
fi

# Mirrors ../data/metadata_8samples.csv exactly: WT infection vs mock, 4h
# and 24h, 2 replicates each.
SRR_LIST=("SRR2936892" "SRR2936893" "SRR2936901" "SRR2936902" "SRR2936907" "SRR2936908" "SRR2936916" "SRR2936917")
GSM_LIST=("GSM1924778" "GSM1924779" "GSM1924787" "GSM1924788" "GSM1924793" "GSM1924794" "GSM1924802" "GSM1924803")
CONDITION_LIST=("infected" "infected" "infected" "infected" "mock" "mock" "mock" "mock")
TIME_LIST=("4h" "4h" "24h" "24h" "4h" "4h" "24h" "24h")

PASSED_SRRS=(); FAILED_SRRS=()

for i in "${!SRR_LIST[@]}"; do
    SRR="${SRR_LIST[$i]}"
    GSM="${GSM_LIST[$i]}"
    CONDITION="${CONDITION_LIST[$i]}"
    TIME="${TIME_LIST[$i]}"
    LOG="logs/${SRR}_pipeline.log"
    echo "Processing ${SRR} (${GSM}, ${CONDITION}, ${TIME})..."

    if ! (
        set -e
        set -o pipefail

        # --- A. Fetch & verify (single-end) ---
        if [ ! -s "raw_reads/${SRR}.fastq" ]; then
            module purge
            module load SRA-Toolkit/3.0.3-gompi-2022a
            prefetch "${SRR}" --output-directory sra_cache >> "${LOG}" 2>&1
            vdb-validate "sra_cache/${SRR}/${SRR}.sra" 2>&1 | grep -q "is consistent"
            fasterq-dump "sra_cache/${SRR}/${SRR}.sra" --outdir raw_reads --threads 8 >> "${LOG}" 2>&1
        fi
        md5sum "raw_reads/${SRR}.fastq" > "checksums/${SRR}_raw.md5"

        # --- B. Raw diagnostic QC ---
        module purge
        module load FastQC/0.11.9-Java-11
        fastqc "raw_reads/${SRR}.fastq" -o qc_before -t 8 --extract >> "${LOG}" 2>&1

        # --- C. Evidence-based fastp trimming (single-end) ---
        module purge
        module load fastp/0.23.4-GCC-13.2.0
        SUMMARY="qc_before/${SRR}_fastqc/summary.txt"
        FASTP_FLAGS=("--length_required" "36")
        if grep -qE "^(WARN|FAIL)\s+Adapter Content" "${SUMMARY}"; then
            FASTP_FLAGS+=("--detect_adapter_for_pe")
        fi
        if grep -qE "^(WARN|FAIL)\s+Per base sequence quality" "${SUMMARY}"; then
            FASTP_FLAGS+=("--cut_tail" "--cut_tail_mean_quality" "20")
        fi
        fastp -i "raw_reads/${SRR}.fastq" -o "trimmed_reads/${SRR}_clean.fastq" \
              "${FASTP_FLAGS[@]}" --thread 8 \
              --json "fastp_reports/${SRR}_fastp.json" --html "fastp_reports/${SRR}_fastp.html" >> "${LOG}" 2>&1
        echo "\`fastp ${FASTP_FLAGS[*]}\`" > "fastp_reports/${SRR}_THRESHOLD_RATIONALE.md"

        # --- D. Dual-reference alignment (single-end) ---
        module purge
        module load STAR/2.7.10b-GCC-11.3.0
        STAR --runThreadN 8 --genomeDir "${REFDIR}/star_index" \
             --readFilesIn "trimmed_reads/${SRR}_clean.fastq" \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix "alignments/${SRR}_" >> "${LOG}" 2>&1

        module purge
        module load SAMtools/1.19.2-GCC-13.2.0
        samtools index "alignments/${SRR}_Aligned.sortedByCoord.out.bam"
        samtools idxstats "alignments/${SRR}_Aligned.sortedByCoord.out.bam" > "qc_reports/${SRR}_idxstats.txt"

        # --- E. Host/pathogen split-quantification ---
        module purge
        module load Subread/2.0.6-GCC-12.3.0
        BAM="alignments/${SRR}_Aligned.sortedByCoord.out.bam"
        featureCounts -T 8 -s 0 -a "${REFDIR}/host_GRCh38.gtf" \
            -o "counts/${SRR}_host_counts.txt" "${BAM}" >> "${LOG}" 2>&1
        featureCounts -T 8 -s 0 -t gene -a "${REFDIR}/pathogen_SL1344.gtf" \
            -o "counts/${SRR}_pathogen_counts.txt" "${BAM}" >> "${LOG}" 2>&1
    ); then
        FAILED_SRRS+=("${SRR}")
        echo "  RESULT: ${SRR} (${GSM}) pipeline FAILED. See ${LOG}."
        continue
    fi

    PASSED_SRRS+=("${SRR}")
    echo "  RESULT: ${SRR} (${GSM}) pipeline PASSED."
done

echo "=========================================================="
echo "Passed: ${PASSED_SRRS[*]:-none}"
echo "Failed: ${FAILED_SRRS[*]:-none}"
echo "=========================================================="

if [ ${#PASSED_SRRS[@]} -eq 0 ]; then
    echo "ERROR: no samples passed. Nothing to merge." >&2
    exit 1
fi

module purge
module load Subread/2.0.6-GCC-12.3.0

HOST_BAMS=()
for SRR in "${PASSED_SRRS[@]}"; do
    HOST_BAMS+=("alignments/${SRR}_Aligned.sortedByCoord.out.bam")
done

featureCounts -T 8 -s 0 -a "${REFDIR}/host_GRCh38.gtf" \
    -o "counts/cohort_host_counts.txt" "${HOST_BAMS[@]}"
featureCounts -T 8 -s 0 -t gene -a "${REFDIR}/pathogen_SL1344.gtf" \
    -o "counts/cohort_pathogen_counts.txt" "${HOST_BAMS[@]}"

echo "OK: counts/cohort_host_counts.txt and counts/cohort_pathogen_counts.txt written."
echo "    These are the two files analysis_dge_host_pathogen/ expects to find."
