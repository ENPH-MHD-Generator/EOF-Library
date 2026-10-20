# Elmer from the base image's source tree, rebuilt with gcc/gfortran 9 against
# MUMPS and Hypre from deps.mk. The build directory is a BuildKit cache mount,
# so rebuilds are incremental.
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

ELMER_DEBUG ?= 0

ELMER_FLAGS := \
  -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran-9 \
  -DCMAKE_C_COMPILER=/usr/bin/gcc-9 \
  -DCMAKE_CXX_COMPILER=/usr/bin/g++-9 \
  -DWITH_MPI=TRUE \
  -DWITH_Mumps=TRUE \
  -DMUMPS_ROOT=$(MUMPS_PREFIX) \
  -DSCALAPACK_LIBRARIES=$(MUMPS_PREFIX)/lib/libscalapack.a \
  -DWITH_Hypre=TRUE \
  -DHYPRE_ROOT=$(HYPRE_PREFIX)

ifeq ($(ELMER_DEBUG),1)
  ELMER_FLAGS += -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow"
else
  ELMER_FLAGS += -DCMAKE_BUILD_TYPE=Release
endif

ELMER_SOLVER_LIB := /usr/local/lib/elmersolver/libelmersolver.so

.PHONY: elmer

elmer:
	export $(GCC9_ENV)
	mkdir -p $(ELMER_SOURCE)/build && cd $(ELMER_SOURCE)/build
	cmake .. $(ELMER_FLAGS)
	# `make install`, not install/fast: the latter installs existing targets
	# without rebuilding, so option changes (MUMPS, debug) would not take effect
	make -j$(JOBS) install
	# Symbols to a file first: under pipefail, grep -q exiting early fails nm
	nm -D $(ELMER_SOLVER_LIB) > /tmp/elmersolver.symbols
	grep -qi " T dmumps" /tmp/elmersolver.symbols \
	  || { echo "Elmer was built without MUMPS; see $(ELMER_SOURCE)/build/CMakeFiles/CMakeError.log" >&2; exit 1; }
	grep -q " T HYPRE_BoomerAMGCreate" /tmp/elmersolver.symbols \
	  || { echo "Elmer was built without Hypre; see $(ELMER_SOURCE)/build/CMakeFiles/CMakeError.log" >&2; exit 1; }
	rm /tmp/elmersolver.symbols
