#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="$SCRIPT_DIR/out"
SAVED_RUNS_DIR="$SCRIPT_DIR/saved_runs"

if [[ $# -ne 1 || "$1" == "-h" || "$1" == "--help" ]]; then
  echo "Usage: $(basename "$0") <name>"
  echo "Saves VTK files from out/ to saved_runs/YYYY-MM-DD_HH-MM_<name>/"
  exit 0
fi

if [[ ! -d "$OUT_DIR" ]]; then
  echo "Error: Output directory not found: $OUT_DIR"
  exit 1
fi

TIMESTAMP=$(date "+%y-%m%d-%H%M")
SAFE_NAME=$(echo "$1" | tr ' ' '_' | tr -cd '[:alnum:]_-')
DEST="$SAVED_RUNS_DIR/${TIMESTAMP}_${SAFE_NAME}"

mkdir -p "$DEST"
find "$OUT_DIR" -name "*.vtk" | while read -r f; do
  GROUP=$(basename "$(dirname "$f")")
  mkdir -p "$DEST/$GROUP"
  cp "$f" "$DEST/$GROUP/"
done

echo "Saved to: $DEST ($(du -sh "$DEST" | cut -f1))"
