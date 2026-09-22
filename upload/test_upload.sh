#!/usr/bin/env bash

FILE="$1"
FILENAME=$(basename "$FILE")

TYPE="local_nexrad_analysis"

RAW_DATETIME="${FILENAME#local_nexrad_analysis_}"
RAW_DATETIME="${RAW_DATETIME%.png}"

DATE="${RAW_DATETIME:0:8}"
TIME="${RAW_DATETIME:9:6}"

DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

echo "Uploading: $FILENAME"
echo "Datetime: $DATETIME"

curl \
  --location \
  "https://weather.atmos.und.edu/api/graphics/upload/${TYPE}/${DATETIME}" \
  --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
  --form "=@${FILE}" \
  --form "fileName=${FILENAME}" \
  --write-out "\nHTTP status: %{http_code}\n"

echo
echo "Upload command finished."