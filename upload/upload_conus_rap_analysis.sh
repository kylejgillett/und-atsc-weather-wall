#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_conus_rap_analysis.py"
STAGE_DIR="staged_figures/conus_rap_analysis"

TYPE="conus_rap_analysis"
BASE_URL="https://weather.atmos.und.edu"

# Clear old CONUS RAP graphics
rm -f "$STAGE_DIR"/conus_analysis_*.png

# Build new CONUS RAP graphics
python "$BUILD_SCRIPT"

# Upload CONUS RAP graphics in map order
for FILE in "$STAGE_DIR"/conus_analysis_*.png; do

    FILENAME=$(basename "$FILE")

    RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}' | head -1)

    DATE="${RAW_DATETIME%_*}"
    TIME="${RAW_DATETIME#*_}"

    DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

    SUFFIX=$(echo "$FILENAME" | sed -E 's/^conus_analysis_[0-9]{8}_[0-9]{6}_([^.]*)\.png$/\1/')

    #echo "Uploading: $FILENAME"
    #echo "Datetime:  $DATETIME"
    #echo "Suffix:    $SUFFIX"

    curl --location \
      "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}/${SUFFIX}" \
      --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
      --form "=@${FILE}" \
      --form "fileName=${FILENAME}" \
      --write-out "\n  + UPLOAD STATUS: http-%{http_code}\n"

done

#echo "CONUS RAP analysis upload finished."