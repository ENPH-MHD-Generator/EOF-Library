# The foundation image built from source (docker/Dockerfile.of6-elmer84): the parts
# of eoflibrary/eof_elmer84_of6 that the solver image uses, on any architecture.
# OpenFOAM 6 is the 20180710 release that image packages, built with Ubuntu
# 16.04's gcc 5 against its Open MPI 1.10.2; Elmer is the same source revision.
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

OPENFOAM_TAG := 20180710
OPENFOAM_URL := https://github.com/OpenFOAM/OpenFOAM-6/archive/refs/tags/$(OPENFOAM_TAG).tar.gz
THIRDPARTY_URL := https://github.com/OpenFOAM/ThirdParty-6/archive/refs/tags/version-6.tar.gz
OPENFOAM_DIR := /opt/openfoam6
THIRDPARTY_DIR := /opt/ThirdParty-6

ELMER_REPO     := https://github.com/ElmerCSC/elmerfem.git
ELMER_REVISION := ed16b68fcb1a761fbe0e28739c5c95c9a34e84aa

.PHONY: packages user openfoam elmer-source

packages:
	apt-get update
	apt-get install -y --no-install-recommends \
	  build-essential gfortran cmake flex libreadline-dev zlib1g-dev \
	  libopenmpi-dev openmpi-bin libblas-dev liblapack-dev libopenblas-dev \
	  git wget ca-certificates python sudo
	rm -rf /var/lib/apt/lists/*

# The base image's user, with the same IDs (user-mapping.sh remaps them at start)
user:
	groupadd -g 98765 openfoam
	useradd -m -u 98765 -g 98765 -G sudo -s /bin/bash openfoam

openfoam:
	mkdir -p $(OPENFOAM_DIR) $(THIRDPARTY_DIR)
	wget -qO- $(OPENFOAM_URL) | tar xz --strip-components=1 -C $(OPENFOAM_DIR)
	wget -qO- $(THIRDPARTY_URL) | tar xz --strip-components=1 -C $(THIRDPARTY_DIR)
	cd $(OPENFOAM_DIR)
	# The package's layout: installed at /opt/openfoam6, user files in $HOME
	sed -i \
	  -e 's|^export FOAM_INST_DIR=$$(cd .*|export FOAM_INST_DIR=/opt|' \
	  -e 's|^export FOAM_INST_DIR=$$HOME/$$WM_PROJECT$$|# &|' \
	  -e 's|^export WM_PROJECT_DIR=.*|export WM_PROJECT_DIR=$$WM_PROJECT_INST_DIR/openfoam6|' \
	  -e 's|^export WM_PROJECT_USER_DIR=.*|export WM_PROJECT_USER_DIR=$$HOME|' \
	  etc/bashrc
	# OpenFOAM 6 predates 64-bit ARM support: add it as OpenFOAM 8 later did,
	# with the x86-64 compiler rules minus -m64
	sed -i 's|^    ppc64le)$$|    aarch64)\n        WM_ARCH=linuxARM64\n        export WM_COMPILER_LIB_ARCH=64\n        export WM_CC=gcc\n        export WM_CXX=g++\n        export WM_CFLAGS=-fPIC\n        export WM_CXXFLAGS="-fPIC -std=c++0x"\n        export WM_LDFLAGS=\n        ;;\n\n&|' \
	  etc/config.sh/settings
	cp -r wmake/rules/linux64Gcc wmake/rules/linuxARM64Gcc
	sed -i 's| -m64||' wmake/rules/linuxARM64Gcc/c wmake/rules/linuxARM64Gcc/c++
	# OpenFOAM's bashrc references unset variables and returns nonzero statuses
	set +eu; source etc/bashrc; set -eu
	test -d wmake/rules/$$WM_ARCH$$WM_COMPILER
	# Scotch and PT-Scotch only; ParaView and CGAL are not used
	cd $(THIRDPARTY_DIR) && ./Allwmake -j $(JOBS) > /tmp/thirdparty.log 2>&1 \
	  || { tail -50 /tmp/thirdparty.log; exit 1; }
	cd $(OPENFOAM_DIR) && ./Allwmake -j $(JOBS) > /tmp/openfoam.log 2>&1 \
	  || { tail -50 /tmp/openfoam.log; exit 1; }
	# Allwmake carries on past failures: check for what the cases use
	for app in rhoPimpleFoam decomposePar reconstructPar foamToVTK potentialFoam gmshToFoam checkMesh; do
	  test -x "$$FOAM_APPBIN/$$app" || { echo "OpenFOAM build is missing $$app" >&2; grep -n -i -B2 -A8 " error" /tmp/openfoam.log | head -80; exit 1; }
	done
	for lib in libfluidThermophysicalModels libcompressibleTurbulenceModels libfiniteVolume libptscotchDecomp libscotchDecomp; do
	  find $$FOAM_LIBBIN -name "$$lib.so" | grep -q . || { echo "OpenFOAM build is missing $$lib" >&2; exit 1; }
	done
	# Object files, dependency files and build trees are not needed at run time
	rm -rf $(OPENFOAM_DIR)/platforms/*/src $(OPENFOAM_DIR)/platforms/*/applications \
	  $(THIRDPARTY_DIR)/build /tmp/thirdparty.log /tmp/openfoam.log

# Only the source tree: elmer.mk configures and builds it with gcc 9
elmer-source:
	mkdir -p $(ELMER_SOURCE) && cd $(ELMER_SOURCE)
	git init -q
	git remote add origin $(ELMER_REPO)
	git fetch -q --depth 1 origin $(ELMER_REVISION)
	git checkout -q FETCH_HEAD
