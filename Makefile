# Host-side build of the EOF-Library MHD solver image. The recipes that run
# inside the image live in docker/mk.
#
#   make solver-image       OpenFOAM + Elmer + EOF coupler + MHD solvers
#   make foundation-image   OpenFOAM 6 and the Elmer source from source, for
#                           architectures without eoflibrary/eof_elmer84_of6
#   make clean              remove the image and its containers
#
# The platform follows the host: linux/arm64 on Apple silicon and other ARM
# hosts, where the solvers then run natively instead of under x86 emulation,
# and linux/amd64 elsewhere. Override with PLATFORM=linux/amd64 (or arm64).

SOLVER_IMAGE     ?= eof-mhd-solvers
FOUNDATION_IMAGE ?= eof-of6-elmer84
ELMER_DEBUG      ?= 0

HOST_ARCH := $(shell uname -m)
ifneq ($(filter arm64 aarch64,$(HOST_ARCH)),)
  PLATFORM ?= linux/arm64
else
  PLATFORM ?= linux/amd64
endif

# eoflibrary/eof_elmer84_of6 is published for linux/amd64 only
ifeq ($(PLATFORM),linux/amd64)
  BASE_IMAGE ?= eoflibrary/eof_elmer84_of6:latest
  BASE_TARGET :=
else
  BASE_IMAGE ?= $(FOUNDATION_IMAGE):latest
  BASE_TARGET := foundation-image
endif

.PHONY: solver-image foundation-image clean

solver-image: $(BASE_TARGET)
	docker build \
	  --platform $(PLATFORM) \
	  --progress=plain \
	  --build-arg BASE_IMAGE=$(BASE_IMAGE) \
	  --build-arg ELMER_DEBUG=$(ELMER_DEBUG) \
	  -f docker/Dockerfile \
	  -t $(SOLVER_IMAGE):latest \
	  .

# Cached by Docker after the first build (about 15 minutes, mostly OpenFOAM)
foundation-image:
	docker build \
	  --platform $(PLATFORM) \
	  --progress=plain \
	  -f docker/Dockerfile.of6-elmer84 \
	  -t $(FOUNDATION_IMAGE):latest \
	  .

clean:
	docker ps -a --filter "ancestor=$(SOLVER_IMAGE)" -q | xargs -r docker rm -f
	docker images $(SOLVER_IMAGE) -q | xargs -r docker rmi -f
