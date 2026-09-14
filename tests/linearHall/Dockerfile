# Linear Hall MHD application: case preparation (caseprep, plasma_collisions,
# Gmsh) and the `mhd` command, on top of the EOF-Library solver image.
# Build context is this directory; build with `make image` here.

ARG SOLVER_IMAGE=eof-mhd-solvers:latest

FROM ghcr.io/astral-sh/uv:0.9.21 AS uv

# Python environment. Ubuntu 16.04 matches the solver image's glibc, so the
# interpreter and wheels resolved here run there.
FROM ubuntu:16.04 AS python
RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates git make \
 && rm -rf /var/lib/apt/lists/*
COPY --from=uv /uv /usr/local/bin/uv
WORKDIR /opt/linear-hall
COPY docker/python-env.mk pyproject.toml uv.lock ./
# plasma_collisions is a private repository: the GitHub token arrives as a build
# secret and is used through a git config in a tmpfs home, never an image layer.
RUN --mount=type=secret,id=github_token,required=true \
    --mount=type=tmpfs,target=/build-home \
    HOME=/build-home make -f python-env.mk python-env


FROM ${SOLVER_IMAGE} AS app
USER root
# Shared-library dependency of the Gmsh wheel
RUN apt-get update \
 && apt-get install -y --no-install-recommends libglu1-mesa \
 && rm -rf /var/lib/apt/lists/*
COPY --from=python /opt/linear-hall/env /opt/linear-hall/env
COPY . /opt/linear-hall
RUN install -m 0755 /opt/linear-hall/mhd-command /usr/local/bin/mhd \
 && install -m 0755 /opt/linear-hall/docker/user-mapping.sh /user-mapping.sh

ENV LINEAR_HALL_HOME=/opt/linear-hall \
    PATH=/opt/linear-hall/env/venv/bin:${PATH} \
    PLASMA_COLLISIONS_BOLSIG_DIR=/opt/linear-hall/env/bolsig \
    PLASMA_COLLISIONS_NO_DOWNLOAD=1 \
    PLASMA_COLLISIONS_CACHE_DIR=/tmp/plasma_collisions

ENTRYPOINT ["/user-mapping.sh"]
