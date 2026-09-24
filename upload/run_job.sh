#!/usr/bin/env bash

# Usage:
#   ./upload/run_job.sh upload_asos.sh
#
#   ./upload/run_job.sh \
#       upload_local_nexrad_analysis.sh \
#       upload_regional_rap_analysis.sh

set -u
set -o pipefail

# Determine repository root automatically from the location of this script.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPLOAD_DIR="${REPO_DIR}/upload"
VENV_DIR="/home/kjgill/weather_wall/undatscweatherwallenv"
ENV_FILE="/home/kjgill/.config/weather-wall/env"
LOG_DIR="/home/kjgill/weather_wall/logs"
MASTER_LOG="${LOG_DIR}/master.log"

source "$ENV_FILE"
export VIRTUAL_ENV="$VENV_DIR"
export PATH="${VENV_DIR}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
cd "$REPO_DIR"

# prevent double or dual running jobs
LOCK_NAME="$(basename "$1" .sh)"
LOCK_FILE="/tmp/weather-wall-${LOCK_NAME}.lock"

exec 9>"$LOCK_FILE"

if ! /usr/bin/flock -n 9; then

    TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

    echo "$TIME_NOW | SKIPPED | $* | already running" >> "$MASTER_LOG"

    echo "$TIME_NOW - Job already running. Skipping."

    exit 0
fi


TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
echo "$TIME_NOW | RUN     | $*" >> "$MASTER_LOG"


STATUS=0
for JOB in "$@"; do
    SCRIPT="${UPLOAD_DIR}/${JOB}"
    echo
    echo "========================================================================"
    echo "------ UND ATSC LEON F OSBORE WEATHER MAPWALL JOB SUBMISSION ------"
    echo "       "
    echo "  + UTC Start: $(date -u '+%Y-%m-%d %H:%M:%SZ')"
    echo "  + Script:    $JOB"
    echo "========================================================================"
    echo

    START_SECONDS=$(date +%s)

    if "$SCRIPT"; then
        END_SECONDS=$(date +%s)
        RUNTIME=$((END_SECONDS - START_SECONDS))

        TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

        echo "$TIME_NOW | SUCCESS | $JOB | ${RUNTIME} sec" \
            >> "$MASTER_LOG"

        echo
        echo "  + Completed: $JOB"
        echo "  + UTC End:   $(date -u '+%Y-%m-%d %H:%M:%SZ')"
    else
        EXIT_CODE=$?

        END_SECONDS=$(date +%s)
        RUNTIME=$((END_SECONDS - START_SECONDS))

        TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

        echo "$TIME_NOW | FAILED  | $JOB | exit=${EXIT_CODE} | ${RUNTIME} sec" \
            >> "$MASTER_LOG"

        echo
        echo "  ! FAILED:    $JOB"
        echo "  ! Exit code: $EXIT_CODE"
        echo "  ! UTC End:   $(date -u '+%Y-%m-%d %H:%M:%SZ')"
        STATUS=$EXIT_CODE
    fi
done


echo
echo "========================================================================"
echo "  + JOB SEQUENCE FINISHED"
echo "  + UTC: $(date -u '+%Y-%m-%d %H:%M:%SZ')"
echo "========================================================================"

exit "$STATUS"