#/bin/bash

DELETE_ORIGINALS=false
if [[ "$1" == "--remove-files" ]]; then
	DELETE_ORIGINALS=true
	echo "Will delete original files."
fi

DEST_DIR="backups"
mkdir -p "$DEST_DIR"

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
FILE_NAME="back-$TIMESTAMP.tar.gz"
TARGETS="tables diffs log latex data"

echo "Compressing into $DEST_DIR/$FILE_NAME..."

if tar -czf "$DEST_DIR/$FILE_NAME" $TARGETS; then
	echo "Success!"

	if [ "$DELETE_ORIGINALS" = true ]; then
		echo "Removing original files..."
		rm -rf $TARGETS
		echo "Done."
	fi
else
	echo "Failed to archive files."
	exit 1
fi
