#!/bin/bash
###############################################################################
# run_pipeline_from_yaml.sh
#
# NARRATIVE:
# Reads pipeline_config.yaml and submits every step it lists as a real
# Slurm job, chained via --dependency=afterok:<previous_jobid>. This is
# Slurm's own native mechanism for "run job B only if job A succeeded" --
# no external workflow-manager software needed, and every step still runs
# with its own script's real #SBATCH resource directives untouched.
#
# Run this ONCE, from a login node -- it returns immediately after
# submitting all jobs. Slurm's scheduler then runs each step in turn,
# entirely unattended, over however many hours/days the real pipeline
# takes.
#
# ISOLATION FROM MANUAL RUNS (RUN_TAG):
# Every numbered script builds its real working directory with an optional
# ${RUN_TAG:+_${RUN_TAG}} suffix -- when RUN_TAG is unset, that expands to
# exactly the original hardcoded path, so a plain `sbatch script.sh` call
# (e.g. a student manually working through the walkthrough) behaves
# identically to before, with zero risk to already-completed manual runs.
# This orchestrator sets a REAL, non-empty RUN_TAG by default
# (`orchestrated`) and exports it into every job's environment via
# `sbatch --export`, so a full automated run always lands in its own
# separate, clean set of directories -- this exact collision (an
# orchestrated run colliding with a completed manual run's output) was
# confirmed as a real failure on a sibling project in this curriculum
# (PRJNA623374: fasterq-dump exited with "file already exists").
# A shared REFDIR (reference genome/index, when this pipeline builds one)
# deliberately does NOT get a RUN_TAG suffix -- it's expensive to rebuild
# and every script already checks for its output before rebuilding it, so
# reusing it across tagged and untagged runs alike is safe.
#
# Usage:
#   bash run_pipeline_from_yaml.sh [path/to/pipeline_config.yaml] \
#        [--start-from N] [--run-tag TAG]
#
#   --start-from N   Skip steps before N and submit only step N onward,
#                     with step N running immediately (no dependency) and
#                     every step after it still chained normally. Use this
#                     to resume after fixing a script that failed partway
#                     through a previous run, without re-submitting steps
#                     that already succeeded.
#
#   --run-tag TAG    Use TAG instead of the default "orchestrated" as the
#                     directory-name suffix for this run. Use a fresh tag
#                     any time you want a genuinely clean slate alongside
#                     a previous orchestrated run's output, without
#                     deleting anything.
#
# WHAT HAPPENS IF A STEP FAILS:
# Slurm's --dependency=afterok means if step N exits non-zero, every
# step after it is automatically left in a "DependencyNeverSatisfied"
# state and never runs -- it does NOT retry, skip, or silently continue.
# Check `sacct` (see this script's own final echo) to see exactly which
# step failed and why before resubmitting.
###############################################################################
set -e
set -o pipefail

YAML_FILE="pipeline_config.yaml"
START_FROM=1
RUN_TAG="orchestrated"

EXPECT_VALUE_FOR=""
for arg in "$@"; do
    if [ -n "${EXPECT_VALUE_FOR}" ]; then
        case "${EXPECT_VALUE_FOR}" in
            start-from) START_FROM="${arg}" ;;
            run-tag)    RUN_TAG="${arg}" ;;
        esac
        EXPECT_VALUE_FOR=""
        continue
    fi
    case "${arg}" in
        --start-from) EXPECT_VALUE_FOR="start-from" ;;
        --run-tag)    EXPECT_VALUE_FOR="run-tag" ;;
        *.yaml|*.yml) YAML_FILE="${arg}" ;;
    esac
done

if [ ! -f "${YAML_FILE}" ]; then
    echo "ERROR: pipeline config not found: ${YAML_FILE}" >&2
    exit 1
fi

SCRIPT_DIR_REL=$(grep -E '^scripts_dir:' "${YAML_FILE}" | sed -E 's/^scripts_dir:\s*//')
YAML_DIR="$(cd "$(dirname "${YAML_FILE}")" && pwd)"
SCRIPTS_DIR="$(cd "${YAML_DIR}/${SCRIPT_DIR_REL}" && pwd)"

echo "Pipeline config   : ${YAML_FILE}"
echo "Scripts dir       : ${SCRIPTS_DIR}"
echo "Starting from step: ${START_FROM}"
echo "Run tag           : ${RUN_TAG} (real working directories get a"
echo "                    _${RUN_TAG} suffix -- isolated from any manual,"
echo "                    untagged run's output)"
echo ""

PREV_JOBID=""
SUBMITTED_ANY=0

while IFS= read -r line; do
    ID=$(echo "${line}" | sed -E 's/.*id:\s*([0-9]+).*/\1/')
    SCRIPT=$(echo "${line}" | sed -E 's/.*script:\s*([^,}]+).*/\1/' | xargs)
    MODE=$(echo "${line}" | sed -E 's/.*mode:\s*([^,}]+).*/\1/' | xargs)

    if [ "${ID}" -lt "${START_FROM}" ]; then
        echo "Skipping step ${ID} (${SCRIPT}) -- before --start-from ${START_FROM}."
        continue
    fi

    SCRIPT_PATH="${SCRIPTS_DIR}/${SCRIPT}"
    if [ ! -f "${SCRIPT_PATH}" ]; then
        echo "ERROR: step ${ID}'s script not found: ${SCRIPT_PATH}" >&2
        exit 1
    fi

    DEP_ARGS=()
    if [ -n "${PREV_JOBID}" ]; then
        DEP_ARGS=(--dependency=afterok:"${PREV_JOBID}")
    fi

    case "${MODE}" in
        sbatch)
            JOBID=$(sbatch --parsable --export=ALL,RUN_TAG="${RUN_TAG}" \
                "${DEP_ARGS[@]}" "${SCRIPT_PATH}")
            ;;
        sbatch_wrap)
            # --partition must be explicit here: --wrap does NOT read the
            # #SBATCH headers from the wrapped script file (there's nothing
            # to read -- the command is a literal string), so nothing
            # supplies a default partition the way it does for a normal
            # `sbatch script.sh` call. Confirmed as a real failure on this
            # cluster ("You must specify a partition (queue) to run in.").
            JOBID=$(sbatch --parsable --export=ALL,RUN_TAG="${RUN_TAG}" \
                "${DEP_ARGS[@]}" \
                --partition=batch --job-name="step${ID}_wrapped" --ntasks=1 \
                --cpus-per-task=1 --mem=2G --time=00:10:00 \
                --wrap="bash ${SCRIPT_PATH}")
            ;;
        *)
            echo "ERROR: unknown mode '${MODE}' for step ${ID}" >&2
            exit 1
            ;;
    esac

    if [ -n "${PREV_JOBID}" ]; then
        echo "Step ${ID} (${SCRIPT}): submitted as job ${JOBID}, depends on job ${PREV_JOBID}"
    else
        echo "Step ${ID} (${SCRIPT}): submitted as job ${JOBID}, no dependency (runs immediately)"
    fi

    PREV_JOBID="${JOBID}"
    SUBMITTED_ANY=1
done < <(grep -E '^\s*-\s*\{' "${YAML_FILE}")

if [ "${SUBMITTED_ANY}" -eq 0 ]; then
    echo "ERROR: no steps submitted -- check --start-from against real step IDs in ${YAML_FILE}." >&2
    exit 1
fi

echo ""
echo "OK: full chain submitted under RUN_TAG=${RUN_TAG}. Monitor with:"
echo "  squeue -u \$(whoami)"
echo "Check final pass/fail per step with:"
echo "  sacct -u \$(whoami) --format=JobID,JobName%20,State,ExitCode,Start,End"
