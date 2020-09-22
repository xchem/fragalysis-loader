#!/usr/bin/env bash

# Does the container look intact?
# We need numerous environment variables.

: "${AWS_ACCESS_KEY_ID?Need to set AWS_ACCESS_KEY_ID}"
: "${AWS_SECRET_ACCESS_KEY?Need to set AWS_SECRET_ACCESS_KEY}"
: "${BUCKET_NAME?Need to set BUCKET_NAME}"
: "${DATA_ORIGIN?Need to set DATA_ORIGIN}"

if [ ! -d "/code/media" ]; then
    echo "ERROR - the destination directory (/code/media) does not exist"
    exit 0
fi

# OK if we get here...
#
# We expect data to be available on S3 in the BUCKET_NAME provided in the
# 'django-data' sub-directory using the path specified in DATA_ORIGIN.
# The bucket data is copied to the mount point '/code/media',
# i.e. data in s3://${BUCKET_NAME)/django-data/${DATA_ORIGIN}
# is written to '/code/media/NEW_DATA'.
#
# - Wipe the (temporary) destination directory
# - Copy new content
# - Run the loader

# Wipe
# ----

DST=/code/media/NEW_DATA
echo "+> Removing ${DST}"
rm -rf ${DST}
mkdir ${DST}

# Copy
# ----

# List the bucket's objects (files).
# Output is typically: -
#
#   2020-09-16 12:56:32      12288 django-data/2020-09-15T16/READY
#   2020-09-16 12:56:32       1886 django-data/2020-09-15T16/TARGET_LIST
#   2020-09-16 12:56:32       4064 django-data/2020-09-15T16/ATAD/ATAD2A-x1712_1/ATAD2A-x1712_1.mol2
#
# And we want...
#
#   django-data/2020-09-15T16/READY
#   django-data/2020-09-15T16/TARGET_LIST
#   django-data/2020-09-15T16/ATAD/ATAD2A-x1712_1/ATAD2A-x1712_1.mol2
BUCKET_PATH="${BUCKET_NAME}/django-data/${DATA_ORIGIN}"
echo "+> Listing S3 path (${BUCKET_PATH})..."
PATH_OBJECTS=$(aws s3 ls --recursive "s3://${BUCKET_PATH}" | tr -s ' ' | cut -d ' ' -f 4)
NUM_PATH_OBJECTS=$(echo $PATH_OBJECTS | tr ' ' '\n' | wc -l)
echo "+> Listed ${NUM_PATH_OBJECTS} objects."

# Now copy each object to the media's NEW_DATA directory.
# The PATH_OBJECT 'django-data/2020-09-15T16/TARGET_LIST'
# is written as 'TARGET_LIST' to '/code/media/NEW_DATA'.
# Do it quietly - we don't want thousands of lines in the log.
echo "+> Copying objects..."
for SRC_OBJECT in $PATH_OBJECTS; do
  DST_OBJECT=$(echo "${SRC_OBJECT}" | cut -f 3- -d '/');
  aws s3 cp "s3://${BUCKET_NAME}/${SRC_OBJECT}" "/code/media/NEW_DATA/${DST_OBJECT}" --only-show-errors;
done
echo "+> Copied."

# Load
# ----

echo "+> Running loader..."
./run_loader.sh
echo "+> Done."
