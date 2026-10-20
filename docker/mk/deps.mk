# Third-party libraries for Elmer, from source with gcc/gfortran 9: ScaLAPACK and
# MUMPS (parallel direct solver) and Hypre. Every Fortran library that shares
# derived types with Elmer must use Elmer's compiler (gfortran 8 changed the
# array descriptor ABI), so Ubuntu's gfortran-5 packages cannot be used.
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

SCALAPACK_VERSION := 2.1.0
SCALAPACK_URL     := https://github.com/Reference-ScaLAPACK/scalapack/archive/refs/tags/v$(SCALAPACK_VERSION).tar.gz
SCALAPACK_SHA256  := f03fda720a152030b582a237f8387014da878b84cbd43c568390e9f05d24617f

MUMPS_VERSION := 5.6.2
MUMPS_URL     := https://mumps-solver.org/MUMPS_$(MUMPS_VERSION).tar.gz
MUMPS_SHA256  := 13a2c1aff2bd1aa92fe84b7b35d88f43434019963ca09ef7e8c90821a8f1d59a

# 2.15.1 matches the age of Elmer 8.4's Hypre interface (Ubuntu's 2.8 is too old)
HYPRE_VERSION := 2.15.1
HYPRE_URL     := https://github.com/hypre-space/hypre/archive/refs/tags/v$(HYPRE_VERSION).tar.gz
HYPRE_SHA256  := 50d0c0c86b4baad227aa9bdfda4297acafc64c3c7256c27351f8bae1ae6f2402

WORK := /tmp/deps-build

# download URL SHA256 FILE
define download
wget -q "$(1)" -O "$(3)"
echo "$(2)  $(3)" | sha256sum -c -
endef

.PHONY: deps scalapack mumps hypre

deps: scalapack mumps hypre

scalapack:
	export $(GCC9_ENV)
	mkdir -p $(WORK) $(MUMPS_PREFIX)/lib && cd $(WORK)
	$(call download,$(SCALAPACK_URL),$(SCALAPACK_SHA256),scalapack.tgz)
	tar xzf scalapack.tgz
	mkdir -p scalapack-$(SCALAPACK_VERSION)/build && cd scalapack-$(SCALAPACK_VERSION)/build
	# The CMake build: the plain Makefile races on the BLACS archive in parallel
	cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_C_COMPILER=mpicc -DCMAKE_Fortran_COMPILER=mpif90 \
	  -DCMAKE_POSITION_INDEPENDENT_CODE=ON -DBUILD_SHARED_LIBS=OFF -DSCALAPACK_BUILD_TESTS=OFF \
	  -DBLAS_LIBRARIES=/usr/lib/libblas.so -DLAPACK_LIBRARIES=/usr/lib/liblapack.so > /dev/null
	make -j$(JOBS) scalapack
	cp lib/libscalapack.a $(MUMPS_PREFIX)/lib/
	rm -rf $(WORK)

mumps: scalapack
	export $(GCC9_ENV)
	mkdir -p $(WORK) $(MUMPS_PREFIX)/lib $(MUMPS_PREFIX)/include && cd $(WORK)
	$(call download,$(MUMPS_URL),$(MUMPS_SHA256),mumps.tgz)
	tar xzf mumps.tgz && cd MUMPS_$(MUMPS_VERSION)
	cp Make.inc/Makefile.inc.generic Makefile.inc
	sed -i \
	  -e 's|^#LMETISDIR .*|LMETISDIR = /usr/lib|' \
	  -e 's|^#IMETIS .*|IMETIS = -I/usr/include|' \
	  -e 's|^#LMETIS .*|LMETIS = -L/usr/lib -lparmetis -lmetis|' \
	  -e 's|^ORDERINGSF .*|ORDERINGSF = -Dpord -Dmetis -Dparmetis|' \
	  -e 's|^CC .*|CC = mpicc|' -e 's|^FC .*|FC = mpif90|' -e 's|^FL .*|FL = mpif90|' \
	  -e 's|^SCALAP .*|SCALAP = $(MUMPS_PREFIX)/lib/libscalapack.a|' \
	  -e 's|^LAPACK .*|LAPACK = -llapack|' -e 's|^LIBBLAS .*|LIBBLAS = -lblas|' \
	  -e 's|^OPTF .*|OPTF = -O3 -fPIC|' -e 's|^OPTC .*|OPTC = -O3 -fPIC -I.|' -e 's|^OPTL .*|OPTL = -O3|' \
	  Makefile.inc
	make -j$(JOBS) d
	cp lib/*.a $(MUMPS_PREFIX)/lib/
	cp include/*.h $(MUMPS_PREFIX)/include/
	rm -rf $(WORK)

# A C library, so no Fortran ABI concerns
hypre:
	export $(GCC9_ENV)
	mkdir -p $(WORK) && cd $(WORK)
	$(call download,$(HYPRE_URL),$(HYPRE_SHA256),hypre.tgz)
	tar xzf hypre.tgz && cd hypre-$(HYPRE_VERSION)/src
	./configure --prefix=$(HYPRE_PREFIX) --with-MPI --disable-fortran CC=mpicc CXX=mpicxx \
	  CFLAGS="-O3 -fPIC" CXXFLAGS="-O3 -fPIC" > /dev/null
	make -j$(JOBS) > /dev/null
	make install > /dev/null
	rm -rf $(WORK)
