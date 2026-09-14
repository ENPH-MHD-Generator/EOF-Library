# Shared settings for the in-image build recipes (docker/Dockerfile runs them).
# Every recipe is one bash script that stops at the first failing command.

SHELL := /bin/bash
.SHELLFLAGS := -euo pipefail -c
.ONESHELL:

# Source trees inside the image
EOF_HOME := /home/openfoam/EOF-Library
EOF_SRC  := $(EOF_HOME)/libs

# Installed locations
OPENFOAM_BASHRC := /opt/openfoam6/etc/bashrc
ELMER_SOURCE    := /opt/elmerfem
# Elmer's cached CMake configuration refers to these prefixes
MUMPS_PREFIX    := /opt/mumps
HYPRE_PREFIX    := /opt/hypre
EOF_LIB         := /opt/eof/lib

# Elmer and its Fortran dependencies use gcc/gfortran 9 (Ubuntu 16.04 ships 5.4);
# OpenFOAM, its coupler and Open MPI keep the system gcc 5. Elmer only uses
# mpif.h, so Open MPI's compiler wrappers are pointed at gcc-9 for those builds.
GCC9_ENV := OMPI_CC=gcc-9 OMPI_CXX=g++-9 OMPI_FC=gfortran-9

JOBS := $(shell nproc)
