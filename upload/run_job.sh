#!/usr/bin/env bash

# Usage:
#   ./upload/run_job.sh upload_asos.sh
#
#   ./upload/run_job.sh \
#       upload_local_nexrad_analysis.sh \
#       upload_regional_rap_analysis.sh

set -u
set -o pipefail


# ============================================================
# PATHS
# ============================================================

# Determine repository root automatically from the location of this script.
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
UPLOAD_DIR="${REPO_DIR}/upload"
VENV_DIR="${HOME}/undatscweatherwallenv"
ENV_FILE="${HOME}/.config/weather-wall/env"
LOG_DIR="${HOME}/weather_wall_logs"
MASTER_LOG="${LOG_DIR}/master.log"

mkdir -p "$LOG_DIR"
echo "$(date -u '+%Y-%m-%dT%H:%M:%SZ') | INVOKED | $*" >> "$MASTER_LOG"


# ============================================================
# ENVIRONMENT
# ============================================================

if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: Environment file not found: $ENV_FILE"
    exit 1
fi

source "$ENV_FILE"

if [ -z "${WEATHER_WALL_API_KEY:-}" ]; then
    echo "ERROR: WEATHER_WALL_API_KEY is not set"
    exit 1
fi

export VIRTUAL_ENV="$VENV_DIR"
export PATH="${VENV_DIR}/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"

cd "$REPO_DIR"


# ============================================================
# HEALTHCHECK
# ============================================================

FIRST_JOB="$(basename "$1")"
case "$FIRST_JOB" in
    upload_local_nexrad_analysis.sh)
        HC_URL="${HC_LOCAL_REGIONAL:-}"
        ;;
    upload_gfs_forecasts.sh)
        HC_URL="${HC_GFS:-}"
        ;;
    upload_soundings_obs.sh)
        HC_URL="${HC_SOUNDINGS_OBS:-}"
        ;;
    upload_conus_rap_analysis.sh)
        HC_URL="${HC_HOURLY_FORECASTS:-}"
        ;;
    upload_asos.sh)
        HC_URL="${HC_ASOS_BUFKIT:-}"
        ;;
    *)
        HC_URL=""
        ;;
esac

# Remove trailing slash if one exists.
HC_URL="${HC_URL%/}"


# ============================================================
# LOCK
# ============================================================

# Prevent double or dual running jobs.
LOCK_NAME="$(basename "$1" .sh)"
LOCK_FILE="/tmp/weather-wall-${LOCK_NAME}.lock"
exec 9>"$LOCK_FILE"
if ! /usr/bin/flock -n 9; then
    TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')
    echo "$TIME_NOW | SKIPPED | $* | already running" >> "$MASTER_LOG"
    echo "$TIME_NOW - Job already running. Skipping."
    exit 0
fi


# ============================================================
# JOB START
# ============================================================

TIME_NOW=$(date -u '+%Y-%m-%dT%H:%M:%SZ')

echo "$TIME_NOW | RUN     | $*" >> "$MASTER_LOG"

# Send Healthchecks START ping.
if [ -n "$HC_URL" ]; then
    curl -fsS \
        -m 10 \
        --retry 3 \
        "${HC_URL}/start" \
        >/dev/null 2>&1 || true
fi


# ============================================================
# RUN JOBS
# ============================================================

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


# ============================================================
# JOB SEQUENCE FINISHED
# ============================================================

echo
echo "========================================================================"
echo "  + JOB SEQUENCE FINISHED"
echo "  + UTC: $(date -u '+%Y-%m-%d %H:%M:%SZ')"
echo "========================================================================"


# ============================================================
# HEALTHCHECK RESULT
# ============================================================

if [ -n "$HC_URL" ]; then
    if [ "$STATUS" -eq 0 ]; then
        # Successful completion.
        curl -fsS \
            -m 10 \
            --retry 3 \
            "$HC_URL" \
            >/dev/null 2>&1 || true
    else
        # One or more scripts in the bundle failed.
        curl -fsS \
            -m 10 \
            --retry 3 \
            "${HC_URL}/fail" \
            >/dev/null 2>&1 || true
    fi
fi


exit "$STATUS"
