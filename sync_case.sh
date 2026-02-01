#!/bin/bash

# Sync case files into a running container (no rebuild).
# Copies case.sif and OpenFOAM case data (0/, constant/, system/, channel.msh)
# into ~/EOF-Library/runs/linearHall so you can edit on the host and run in the container.
# constant/polyMesh is never copied — blockMesh writes it in the container; we only sync g, transportProperties, turbulenceProperties.

set -e

# Find the running eof_local container
CONTAINER_ID=$(docker ps --filter "ancestor=eof_local:latest" --format "{{.ID}}" | head -n 1)

if [ -z "$CONTAINER_ID" ]; then
    echo "Error: No running container found from image eof_local:latest"
    echo "Please start the container first with ./launch.sh"
    exit 1
fi

echo "Found running container: $CONTAINER_ID"

# Project root and case directory (default: linearHall)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CASE_NAME="${1:-linearHall}"
CASE_DIR="$SCRIPT_DIR/tests/$CASE_NAME"
RUN_DIR="/home/openfoam/EOF-Library/runs/$CASE_NAME"

if [ ! -d "$CASE_DIR" ]; then
    echo "Error: Case directory not found: $CASE_DIR"
    exit 1
fi

if [ ! -f "$CASE_DIR/case.sif" ]; then
    echo "Error: case.sif not found in $CASE_DIR"
    exit 1
fi

echo "Syncing case '$CASE_NAME' from $CASE_DIR to container $RUN_DIR ..."

# Ensure run directory exists in container
docker exec "$CONTAINER_ID" mkdir -p "$RUN_DIR"

# Elmer input
docker cp "$CASE_DIR/case.sif" "$CONTAINER_ID:$RUN_DIR/case.sif"

# OpenFOAM initial conditions and config (no rebuild needed)
[ -d "$CASE_DIR/0" ]       && docker cp "$CASE_DIR/0"       "$CONTAINER_ID:$RUN_DIR/"
# constant/ but skip polyMesh — blockMesh writes constant/polyMesh/ in the container; we never overwrite it
if [ -d "$CASE_DIR/constant" ]; then
  docker exec "$CONTAINER_ID" mkdir -p "$RUN_DIR/constant"
  for item in "$CASE_DIR/constant"/*; do
    [ -e "$item" ] || continue
    [ "$(basename "$item")" = "polyMesh" ] && continue
    docker cp "$item" "$CONTAINER_ID:$RUN_DIR/constant/"
  done
fi
[ -d "$CASE_DIR/system" ]   && docker cp "$CASE_DIR/system"   "$CONTAINER_ID:$RUN_DIR/"

# Mesh input for ElmerGrid (needed if you run generate_mesh again)
[ -f "$CASE_DIR/channel.msh" ] && docker cp "$CASE_DIR/channel.msh" "$CONTAINER_ID:$RUN_DIR/channel.msh"

echo ""
echo "Case synced successfully!"
echo "In the container, run from: $RUN_DIR"
echo "  configure_elmer   # if you need to set ELMERSOLVER_STARTINFO"
echo "  run_sim           # to run the coupled simulation"
