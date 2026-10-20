#!/usr/bin/env bash

set -Eeuo pipefail

CONTAINER_USER="openfoam"
CONTAINER_HOME="/home/openfoam"

if [[ $# -eq 0 ]]; then
  set -- /bin/bash
fi

if [[ -n "${HOST_USER_ID:-}" || -n "${HOST_USER_GID:-}" ]]; then
  if [[ -z "${HOST_USER_ID:-}" || -z "${HOST_USER_GID:-}" ]]; then
    echo "HOST_USER_ID and HOST_USER_GID must be supplied together." >&2
    exit 2
  fi
  if [[ ! "$HOST_USER_ID" =~ ^[0-9]+$ || ! "$HOST_USER_GID" =~ ^[0-9]+$ ]]; then
    echo "HOST_USER_ID and HOST_USER_GID must be numeric." >&2
    exit 2
  fi

  current_uid="$(id -u "$CONTAINER_USER")"
  current_gid="$(id -g "$CONTAINER_USER")"
  if [[ "$current_gid" != "$HOST_USER_GID" ]]; then
    groupmod --non-unique --gid "$HOST_USER_GID" "$CONTAINER_USER"
  fi
  if [[ "$current_uid" != "$HOST_USER_ID" || "$current_gid" != "$HOST_USER_GID" ]]; then
    usermod --non-unique --uid "$HOST_USER_ID" --gid "$HOST_USER_GID" "$CONTAINER_USER"
    chown -R "$HOST_USER_ID:$HOST_USER_GID" "$CONTAINER_HOME"
  fi
else
  echo "Warning: HOST_USER_ID and HOST_USER_GID are unset; /runs may contain container-owned files." >&2
fi

exec /sbin/runuser -u "$CONTAINER_USER" -- "$@"
