#!/usr/bin/env bash
set -e

BASE_URL="https://weather.atmos.und.edu"
ENV_FILE="${HOME}/.config/weather-wall/env"

if [ "$#" -ne 2 ]; then
    echo "Usage:"
    echo "  $0 <graphic_type> <earliest_datetime_to_keep>"
    echo
    echo "Example:"
    echo "  $0 local_nexrad_analysis 2026-09-30T18:00:00Z"
    exit 1
fi

TYPE="$1"
DATETIME="$2"

if [ ! -f "$ENV_FILE" ]; then
    echo "ERROR: Environment file not found: $ENV_FILE"
    exit 1
fi

source "$ENV_FILE"

if [ -z "${WEATHER_WALL_API_KEY:-}" ]; then
    echo "ERROR: WEATHER_WALL_API_KEY is not set"
    exit 1
fi

echo
echo "Graphic type:       $TYPE"
echo "Earliest to keep:   $DATETIME"
echo
echo "This will DELETE all '$TYPE' images older than:"
echo "  $DATETIME"
echo

read -r -p "Continue? [y/N] " CONFIRM

if [[ ! "$CONFIRM" =~ ^[Yy]$ ]]; then
    echo "Cancelled."
    exit 0
fi

curl --fail-with-body \
    --location \
    --request DELETE \
    "${BASE_URL}/api/graphics/${TYPE}/${DATETIME}" \
    --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
    --write-out "\nHTTP status: %{http_code}\n"

echo
echo "Delete request completed."