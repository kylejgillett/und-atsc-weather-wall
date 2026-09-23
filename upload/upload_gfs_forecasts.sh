#!/usr/bin/env bash

set -e
shopt -s nullglob

BUILD_SCRIPT="build_figures/build_conus_gfs_forecast.py"
STAGE_DIR="staged_figures/conus_gfs_forecasts"

BASE_URL="https://weather.atmos.und.edu"

# Clear old GFS graphics
rm -f "$STAGE_DIR"/gfs_*.png

# Build new GFS graphics
echo "Building GFS graphics..."
python "$BUILD_SCRIPT"
echo "GFS build finished."

# Upload all GFS graphics
for FILE in "$STAGE_DIR"/gfs_*.png; do

    FILENAME=$(basename "$FILE")

    # Extract type, e.g. gfs_000a, gfs_300a, gfs_500b
    TYPE=$(echo "$FILENAME" | sed -E \
      's/^(gfs_[0-9]{3}[ab])_.*/\1/')

    # Extract cycle datetime
    RAW_DATETIME=$(echo "$FILENAME" | grep -oE \
      '[0-9]{8}_[0-9]{6}' | head -1)

    DATE="${RAW_DATETIME%_*}"
    TIME="${RAW_DATETIME#*_}"

    DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

    # Forecast hour is the suffix
    SUFFIX=$(echo "$FILENAME" | sed -E \
      's/^gfs_[0-9]{3}[ab]_[0-9]{8}_[0-9]{6}_([^.]*)\.png$/\1/')

    echo
    echo "Uploading: $FILENAME"
    echo "Type:      $TYPE"
    echo "Datetime:  $DATETIME"
    echo "Suffix:    $SUFFIX"

    curl --location \
      "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}/${SUFFIX}" \
      --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
      --form "=@${FILE}" \
      --form "fileName=${FILENAME}" \
      --write-out "\nHTTP status: %{http_code}\n"

done

echo
echo "GFS upload finished."
