#!/bin/bash

# Default debug flag
ELMER_DEBUG=0

# Help
usage() {
  cat <<EOF
Usage: $(basename "$0") [options]

Options:
  -d, --debug        Build Elmer with debug flags
  -h, --help         Show this help message and exit

Examples:
  $(basename "$0")
  $(basename "$0") --debug
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
  case "$1" in
    -d|--debug)
      ELMER_DEBUG=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1"
      echo
      usage
      exit 1
      ;;
  esac
done

# Go to the directory that this script is in
cd "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")"

# Build, and exit strictly on failure
make build ELMER_DEBUG=${ELMER_DEBUG} || exit 1

docker run --rm -it \
  --platform linux/amd64 \
  -e HOST_USER_ID=$(id -u) \
  -e HOST_USER_GID=$(id -g) \
  -v "$(pwd)/out:/home/openfoam/EOF-Library/runs" \
  eof_local:latest

# docker run --rm -it --platform linux/amd64 -e HOST_USER_ID=$(id -u) -e HOST_USER_GID=$(id -g) eof_local:latest

