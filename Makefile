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

eof: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean $(EOF_SRC)/coupleElmer
	. $(OPENFOAM_HOME)/etc/bashrc && wmake $(EOF_SRC)/coupleElmer
	elmerf90 -o $(EOF_SRC)/Elmer2OpenFOAM.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/Elmer2OpenFOAM.F90
	elmerf90 -o $(EOF_SRC)/OpenFOAM2Elmer.so -J $(nproc) $(EOF_SRC) $(EOF_SRC)/OpenFOAM2Elmer.F90

solver: environment
	. $(OPENFOAM_HOME)/etc/bashrc && wclean solvers/mdhLinearHall
	. $(OPENFOAM_HOME)/etc/bashrc && wmake solvers/mdhLinearHall
	rm -rf solvers/mdhLinearHall/processor*

# Elmer debug flag
ifeq ($(ELMER_DEBUG),1)
  ELMER_CMAKE_FLAGS := \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow"
else
  ELMER_CMAKE_FLAGS := -DCMAKE_BUILD_TYPE=Release
endif

elmer: environment
	sudo rm -f /opt/elmerfem/fem/src/modules/MHDSolve.F90
	sudo mkdir -p /opt/elmerfem/fem/src/modules/MHDSolve
	sudo cp libs/solvers/MHDSolve/MHDUtils.F90 /opt/elmerfem/fem/src/modules/MHDSolve/
	sudo cp libs/solvers/MHDSolve/MHDSolve.F90 /opt/elmerfem/fem/src/modules/MHDSolve/
	cd /opt/elmerfem/build && sudo cmake .. $(ELMER_CMAKE_FLAGS)
	cd /opt/elmerfem/build && sudo make -j"$$(nproc)" MHDSolve
	cd /opt/elmerfem/build && sudo make install/fast
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
