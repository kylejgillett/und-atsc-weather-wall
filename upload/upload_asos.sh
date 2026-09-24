#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_asos_timeseries.py"
STAGE_DIR="staged_figures/asos_timeseries"

TYPE="asos_timeseries"
BASE_URL="https://weather.atmos.und.edu"

# Clear old ASOS graphics
rm -f "$STAGE_DIR"/asos_timeseries_*.png

# Build new ASOS graphics
python "$BUILD_SCRIPT"

# Upload ASOS graphics in loop order
for FILE in "$STAGE_DIR"/asos_timeseries_*_[0-9][0-9]-*.png; do

    FILENAME=$(basename "$FILE")

    RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}' | head -1)

    DATE="${RAW_DATETIME%_*}"
    TIME="${RAW_DATETIME#*_}"

    DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

    SUFFIX=$(echo "$FILENAME" | sed -E 's/^asos_timeseries_[0-9]{8}_[0-9]{6}_([^.]*)\.png$/\1/')

    # echo "Uploading: $FILENAME"
    # echo "Suffix:    $SUFFIX"

    curl --location \
      "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}/${SUFFIX}" \
      --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
      --form "=@${FILE}" \
      --form "fileName=${FILENAME}" \
      --write-out "\n  + UPLOAD STATUS: http-%{http_code}\n"

done

# echo "ASOS upload finished."