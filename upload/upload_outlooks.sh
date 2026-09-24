#!/usr/bin/env bash

set -e

BUILD_SCRIPT="build_figures/build_outlooks.py"
STAGE_DIR="staged_figures/conus_outlooks"

TYPE="outlook"
BASE_URL="https://weather.atmos.und.edu"

# Clear old outlook graphics
rm -f "$STAGE_DIR"/outlook_*.png

# Build new outlook graphics
python "$BUILD_SCRIPT"

# Upload outlook graphics
for FILE in "$STAGE_DIR"/outlook_*_[0-9][0-9]-*.png; do

    # Skip if no matching files were generated
    [ -e "$FILE" ] || continue

    FILENAME=$(basename "$FILE")

    RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}' | head -1)

    DATE="${RAW_DATETIME%_*}"
    TIME="${RAW_DATETIME#*_}"

    DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

    SUFFIX=$(echo "$FILENAME" | sed -E 's/^outlook_[0-9]{8}_[0-9]{6}_([0-9]{2}-.+)\.png$/\1/')

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

#echo "Outlook upload finished."






























# #!/usr/bin/env bash

# set -e

# BUILD_SCRIPT="build_figures/build_conus_spc_outlooks.py"
# STAGE_DIR="staged_figures/conus_spc_outlooks"

# TYPE="outlook"
# BASE_URL="https://weather.atmos.und.edu"

# # Clear old outlook graphics
# rm -f "$STAGE_DIR"/outlook_*.png

# # Build new outlook graphics
# python "$BUILD_SCRIPT"

# # Upload outlook graphics in loop order
# for FILE in "$STAGE_DIR"/outlook_*_[0-9][0-9].png; do

#     FILENAME=$(basename "$FILE")

#     RAW_DATETIME=$(echo "$FILENAME" | grep -oE '[0-9]{8}_[0-9]{6}' | head -1)

#     DATE="${RAW_DATETIME%_*}"
#     TIME="${RAW_DATETIME#*_}"

#     DATETIME="${DATE:0:4}-${DATE:4:2}-${DATE:6:2}T${TIME:0:2}:${TIME:2:2}:${TIME:4:2}Z"

#     SUFFIX=$(echo "$FILENAME" | sed -E 's/^outlook_[0-9]{8}_[0-9]{6}_([0-9]{2})\.png$/\1/')

#     echo "Uploading: $FILENAME"
#     echo "Datetime:  $DATETIME"
#     echo "Suffix:    $SUFFIX"

#     curl --location \
#       "${BASE_URL}/api/graphics/upload/${TYPE}/${DATETIME}/${SUFFIX}" \
#       --header "X-API-Key: ${WEATHER_WALL_API_KEY}" \
#       --form "=@${FILE}" \
#       --form "fileName=${FILENAME}" \
#       --write-out "\nHTTP status: %{http_code}\n"

# done

# echo "SPC outlook upload finished."