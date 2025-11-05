#!/bin/bash
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <list_file>"
    exit 1
fi

LIST_FILE="$1"

if [ ! -f "$LIST_FILE" ]; then
    echo "Error: list file $LIST_FILE not found"
    exit 1
fi

echo "Start deleting directories from $LIST_FILE..."

while IFS= read -r dir; do
    # 跳过空行
    [ -z "$dir" ] && continue

    if [ -d "$dir" ]; then
        echo "Deleting $dir"
        rm -rf "$dir"
    else
        echo "Skip (not found): $dir"
    fi
done < "$LIST_FILE"

echo "Finished $LIST_FILE"