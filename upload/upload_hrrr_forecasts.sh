#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_hrrr_forecast.py"
STAGE_DIR="staged_figures/hrrr_forecasts"

TYPE="hrrr_forecast"
BASE_URL="https://weather.atmos.und.edu"

rm -f "$STAGE_DIR"/hrrr_forecast_*.png

#Build new HRRR forecast
python "$BUILD_SCRIPT"

# Find HRRR forecast files
FILES=$(ls -1 "$STAGE_DIR"/hrrr_forecast_*.png | sort)

for FILE in $FILES; do

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

done

echo "HRRR forecast upload finished."