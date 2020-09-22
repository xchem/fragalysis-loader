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
# On S3 the object is always called media.tar.gz
# so we pull this file down and decompress it into '/code/media/NEW_DATA'
# before deleting the downloaded media archive.
#
# And, when copied into the stack image, the directories resemble...
#
#   /code/media/NEW_DATA/READY
#   /code/media/NEW_DATA/TARGET_LIST
#   /code/media/NEW_DATA/ATAD/ATAD2A-x1712_1/ATAD2A-x1712_1.mol2
#
BUCKET_PATH="${BUCKET_NAME}/django-data/${DATA_ORIGIN}/media.tar.gz"
echo "+> Getting S3 object (${BUCKET_PATH})..."
aws s3 cp "s3://${BUCKET_PATH}" "/code/media/NEW_DATA/media.tar.gz" --only-show-errors
echo "+> Unpacking..."
tar -zxf /code/media/NEW_DATA/media.tar.gz -C /code/media/NEW_DATA
rm /code/media/NEW_DATA/media.tar.gz
echo "+> Unpacked."

# Load
# ----

echo "+> Running loader..."
./run_loader.sh
echo "+> Done."
