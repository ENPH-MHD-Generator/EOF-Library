export ELMER_HOME        := /usr/local
export ELMER_SOLVER_HOME := /usr/local
export EOF_HOME          := /home/openfoam/EOF-Library
export EOF_SRC           := $(EOF_HOME)/libs
export PATH              := /usr/local/bin:$(PATH)
export LD_LIBRARY_PATH   := /usr/local/lib:$(LD_LIBRARY_PATH)
export OPENFOAM_HOME	 := /opt/openfoam6
export ROOT_DIR			 := $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST)))) # https://stackoverflow.com/a/23324703
export IMAGE_NAME		 ?= mhd-sim

# ---- Build options ----
ELMER_DEBUG ?= 0

# This is so that the environment variables persist between commands
SHELL := /bin/bash
.ONESHELL:

# -- Docker Container

environment:
	. $(OPENFOAM_HOME)/etc/bashrc
	. $(EOF_HOME)/etc/bashrc
	cd $(EOF_HOME)

eof: eof-openfoam eof-elmer

# OpenFOAM side of the coupler (C++, independent of the Elmer build)
eof-openfoam: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean $(EOF_SRC)/coupleElmer
	. $(OPENFOAM_HOME)/etc/bashrc && wmake $(EOF_SRC)/coupleElmer

# Elmer-side Fortran modules. They must be compiled against the installed
# Elmer: its derived types change with build options (e.g. HAVE_MUMPS adds
# matrix fields), so modules built against another Elmer corrupt memory.
eof-elmer: environment
	set -e
	elmerf90 -o $(EOF_SRC)/Elmer2OpenFOAM.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/Elmer2OpenFOAM.F90
	elmerf90 -o $(EOF_SRC)/OpenFOAM2Elmer.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/OpenFOAM2Elmer.F90
	elmerf90 -o $(EOF_SRC)/MHDSolve.so       -J $(nproc) $(EOF_SRC) $(EOF_SRC)/solvers/MHDSolve/MHDUtils.F90 $(EOF_SRC)/solvers/MHDSolve/MHDSolve.F90

solver: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean solvers/mdhLinearHall
	. $(OPENFOAM_HOME)/etc/bashrc && wmake solvers/mdhLinearHall
	rm -rf solvers/mdhLinearHall/processor*

# Elmer with gcc/gfortran 9 and the MUMPS parallel direct solver built from source
# in /opt/mumps (see docker/Dockerfile.build_simulation)
ELMER_SOLVER_FLAGS := \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-9 \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-9 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-9 \
  -DWITH_MPI=TRUE \
  -DWITH_Mumps=TRUE \
  -DMUMPS_ROOT=/opt/mumps \
  -DSCALAPACK_LIBRARIES=/opt/mumps/lib/libscalapack.a \
  -DWITH_Hypre=TRUE \
  -DHYPRE_ROOT=/opt/hypre
ELMER_BUILD_ENV := OMPI_CC=gcc-9 OMPI_CXX=g++-9 OMPI_FC=gfortran-9

# Elmer debug flag
ifeq ($(ELMER_DEBUG),1)
  ELMER_CMAKE_FLAGS := \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow"
else
  ELMER_CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release
endif

elmer: environment
	set -e
	cd /opt/elmerfem/build && sudo env $(ELMER_BUILD_ENV) cmake .. $(ELMER_CMAKE_FLAGS) $(ELMER_SOLVER_FLAGS)
	# install/fast only installs existing targets without building, so flag
	# changes (MUMPS, debug) would silently not take effect
	cd /opt/elmerfem/build && sudo env $(ELMER_BUILD_ENV) make -j$$(nproc) install
	nm -D /usr/local/lib/elmersolver/libelmersolver.so | grep -qi " T dmumps" \
	  || { echo "Elmer was built without MUMPS; see /opt/elmerfem/build/CMakeFiles/CMakeError.log" >&2; exit 1; }
	nm -D /usr/local/lib/elmersolver/libelmersolver.so | grep -q " T HYPRE_BoomerAMGCreate" \
	  || { echo "Elmer was built without Hypre; see /opt/elmerfem/build/CMakeFiles/CMakeError.log" >&2; exit 1; }
	cd $(EOF_HOME)

# -- Host System

build_environment:
	cd $(ROOT_DIR)

setup: build_environment
	mkdir -p ./experiments ./out

build: setup
	docker build \
	  --build-arg ELMER_DEBUG=$(ELMER_DEBUG) \
	  --progress=plain \
	  --network host \
	  --platform linux/amd64 \
	  -f docker/Dockerfile.build_simulation \
	  -t $(IMAGE_NAME):latest .

clean: build_environment
	docker ps -a --filter "ancestor=$(IMAGE_NAME)" -q | xargs -r docker rm -f
	docker images $(IMAGE_NAME) -q | xargs -r docker rmi -f
