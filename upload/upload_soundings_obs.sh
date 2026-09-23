#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_obs_soundings.py"
STAGE_DIR="staged_figures/soundings"

TYPE="sounding"
BASE_URL="https://weather.atmos.und.edu"

# Clear old observed sounding graphics
rm -f "$STAGE_DIR"/sounding_*_[0-9][0-9]-obs-*.png

# Build new observed sounding graphics
python "$BUILD_SCRIPT"

# Upload observed soundings in loop order
for FILE in "$STAGE_DIR"/sounding_*_[0-9][0-9]-obs-*.png; do

    FILENAME=$(basename "$FILE")

    RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}' | head -1)

    DATE="${RAW_DATETIME%_*}"
    TIME="${RAW_DATETIME#*_}"

    DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

    SUFFIX=$(echo "$FILENAME" | sed -E 's/^sounding_[0-9]{8}_[0-9]{6}_([^.]*)\.png$/\1/')

    echo "Uploading: $FILENAME"
    echo "Datetime:  $DATETIME"
    echo "Suffix:    $SUFFIX"

    curl --location \
      "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}/${SUFFIX}" \
      --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
      --form "=@${FILE}" \
      --form "fileName=${FILENAME}" \
      --write-out "\nHTTP status: %{http_code}\n"

done

echo "Observed sounding upload finished."