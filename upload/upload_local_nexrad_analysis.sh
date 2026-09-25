#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_local_nexrad_analysis.py"
STAGE_DIR="staged_figures/local_nexrad_analysis"

TYPE="local_nexrad_analysis"
BASE_URL="https://weather.atmos.und.edu"

rm -f "$STAGE_DIR"/local_nexrad_analysis_*.png

python "$BUILD_SCRIPT"


# Find newest generated image
FILE=$(ls -1 "$STAGE_DIR"/local_nexrad_analysis_*.png | sort | tail -1)
FILENAME=$(basename "$FILE")

# Extract datetime from filename
RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}')

DATE="${RAW_DATETIME%_*}"
TIME="${RAW_DATETIME#*_}"

DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

#echo
#echo "Uploading: $FILENAME"
#echo "Datetime:  $DATETIME"

curl --location \
  "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}" \
  --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
  --form "=@${FILE}" \
  --form "fileName=${FILENAME}" \
  --write-out "\n  + UPLOAD STATUS: http-%{http_code}\n"

#echo
#echo "Local NEXRAD upload finished."
