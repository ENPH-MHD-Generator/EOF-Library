#!/bin/bash

# Quick rebuild script for MHDSolve.F90
# Copies the file to a running container and rebuilds only that module

set -e

# Find the running eof_local container
CONTAINER_ID=$(docker ps --filter "ancestor=eof_local:latest" --format "{{.ID}}" | head -n 1)

if [ -z "$CONTAINER_ID" ]; then
    echo "Error: No running container found from image eof_local:latest"
    echo "Please start the container first with ./launch.sh"
    exit 1
fi

echo "Found running container: $CONTAINER_ID"

# Get the project root directory (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MHDSOLVE_DIR="$SCRIPT_DIR/libs/solvers/MHDSolve"

if [ ! -f "$MHDSOLVE_DIR/MHDSolve.F90" ] || [ ! -f "$MHDSOLVE_DIR/MHDUtils.F90" ]; then
    echo "Error: MHDSolve.F90 and/or MHDUtils.F90 not found in $MHDSOLVE_DIR"
    exit 1
fi

echo "Copying MHDSolve sources to container..."
docker cp "$MHDSOLVE_DIR/MHDUtils.F90" "$CONTAINER_ID:/home/openfoam/EOF-Library/libs/solvers/MHDSolve/MHDUtils.F90"
docker cp "$MHDSOLVE_DIR/MHDSolve.F90" "$CONTAINER_ID:/home/openfoam/EOF-Library/libs/solvers/MHDSolve/MHDSolve.F90"

echo "Rebuilding MHDSolve in container..."
docker exec -it "$CONTAINER_ID" bash -c "
    source /home/openfoam/.bashrc && \
    sudo rm -f /opt/elmerfem/fem/src/modules/MHDSolve.F90 && \
    sudo mkdir -p /opt/elmerfem/fem/src/modules/MHDSolve && \
    sudo cp ~/EOF-Library/libs/solvers/MHDSolve/MHDUtils.F90 /opt/elmerfem/fem/src/modules/MHDSolve/ && \
    sudo cp ~/EOF-Library/libs/solvers/MHDSolve/MHDSolve.F90 /opt/elmerfem/fem/src/modules/MHDSolve/ && \
    cd /opt/elmerfem/build && \
    sudo cmake .. && \
    sudo make MHDSolve && \
    sudo make install MHDSolve
"

echo ""
echo "MHDSolve rebuilt successfully!"
echo "The updated solver is now available in the running container."
