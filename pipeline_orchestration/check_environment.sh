#!/bin/bash
###############################################################################
# check_environment.sh
#
# NARRATIVE:
# A real orchestrated run on a sibling repo in this curriculum
# (PRJNA274277) failed partway through a chain because a script's
# hardcoded module version didn't exist on this cluster's module tree --
# discovered only after earlier steps had already burned real job time.
# This script checks for that whole category of failure BEFORE submitting
# the real chain:
#   1) every `module load` referenced in ../pipeline_raw_reads_to_counts/
#      (the manual walkthrough) and ./scripts/ (this orchestrator's own
#      copies) is checked against what this cluster's module tree actually
#      has (via `module spider`, which resolves modules without loading
#      them or triggering a GCC toolchain swap the way `module load` does);
#   2) the manual walkthrough scripts and the orchestrator's own copies are
#      compared to make sure they still load the SAME modules -- these are
#      deliberately kept as separate files (see this folder's own
#      README.md), so a fix applied to one and forgotten in the other is a
#      real risk this check is meant to catch;
#   3) every reference download URL (FASTA/GTF) referenced by either copy
#      is checked for reachability.
#
# This does NOT try to auto-fix anything, and it will not catch every
# possible failure (e.g. a real logic bug in a script) -- it only catches
# the "the environment doesn't have what this script assumes it has"
# category.
#
# Usage:
#   bash check_environment.sh
#
# Exits 0 if everything checked out, 1 if anything is missing, unreachable,
# or inconsistent -- safe to use as a gate before the real chain:
#   bash check_environment.sh && bash run_pipeline_from_yaml.sh pipeline_config.yaml
###############################################################################
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MANUAL_DIR="${SCRIPT_DIR}/../pipeline_raw_reads_to_counts"
ORCH_DIR="${SCRIPT_DIR}/scripts"

FAIL=0

echo "=== Checking module consistency between manual and orchestrator script copies ==="
for manual_script in "${MANUAL_DIR}"/*.sh; do
    base=$(basename "${manual_script}")
    orch_script="${ORCH_DIR}/${base}"
    [ -f "${orch_script}" ] || continue

    manual_mods=$(grep -oE 'module load [^ ]+' "${manual_script}" | sort -u)
    orch_mods=$(grep -oE 'module load [^ ]+' "${orch_script}" | sort -u)

    if [ "${manual_mods}" != "${orch_mods}" ]; then
        echo "DRIFT: ${base} loads different modules in manual vs orchestrator copy:"
        diff <(echo "${manual_mods}") <(echo "${orch_mods}") | sed 's/^/    /'
        FAIL=1
    fi
done
[ "${FAIL}" -eq 0 ] && echo "OK: no drift detected between manual and orchestrator script copies."
echo ""

echo "=== Checking every referenced module exists on this cluster ==="
ALL_MODULES=$(grep -hoE 'module load [^ ]+' "${MANUAL_DIR}"/*.sh "${ORCH_DIR}"/*.sh | sed 's/module load //' | sort -u)
for mod in ${ALL_MODULES}; do
    OUTPUT=$(module spider "${mod}" 2>&1)
    if echo "${OUTPUT}" | grep -qi "error\|unknown"; then
        echo "MISSING: ${mod}"
        echo "${OUTPUT}" | sed 's/^/    /'
        FAIL=1
    else
        echo "OK: ${mod}"
    fi
done
echo ""

echo "=== Checking reference download URLs are reachable ==="
ALL_URLS=$(grep -hoE '[A-Z_]*URL="[^"]+"' "${MANUAL_DIR}"/*.sh "${ORCH_DIR}"/*.sh | sed -E 's/^[A-Z_]+="//; s/"$//' | sort -u)
if [ -z "${ALL_URLS}" ]; then
    echo "(none referenced by this pipeline)"
else
    for url in ${ALL_URLS}; do
        if wget --spider -q --timeout=20 "${url}"; then
            echo "OK: ${url}"
        else
            echo "UNREACHABLE: ${url}"
            FAIL=1
        fi
    done
fi
echo ""

if [ "${FAIL}" -eq 0 ]; then
    echo "=== ALL CHECKS PASSED -- safe to run: bash run_pipeline_from_yaml.sh pipeline_config.yaml ==="
    exit 0
else
    echo "=== ONE OR MORE CHECKS FAILED -- fix the items above before running the real chain. ===" >&2
    exit 1
fi
