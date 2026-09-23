#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_regional_sfc_analysis.py"
STAGE_DIR="staged_figures/regional_surface_analysis"

TYPE="regional_surface_analysis"
BASE_URL="https://weather.atmos.und.edu"

rm -f "$STAGE_DIR"/regional_rap_analysis_*.png

#Build new RAP analysis
python "$BUILD_SCRIPT"

# Find newest RAP analysis
FILE=$(ls -1 "$STAGE_DIR"/regional_rap_analysis_*.png | sort | tail -1)
FILENAME=$(basename "$FILE")

RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}')

DATE="${RAW_DATETIME%_*}"
TIME="${RAW_DATETIME#*_}"

DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

echo "Uploading: $FILENAME"
echo "Datetime:  $DATETIME"

curl --location \
  "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}" \
  --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
  --form "=@${FILE}" \
  --form "fileName=${FILENAME}" \
  --write-out "\nHTTP status: %{http_code}\n"

echo "Regional RAP analysis upload finished."