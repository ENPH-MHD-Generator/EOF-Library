# Host-side build of the EOF-Library MHD solver image. The recipes that run
# inside the image live in docker/mk.
#
#   make solver-image   OpenFOAM + Elmer + EOF coupler + MHD solvers
#   make clean          remove the image and its containers

SOLVER_IMAGE ?= eof-mhd-solvers
ELMER_DEBUG  ?= 0
PLATFORM     := linux/amd64

.PHONY: solver-image clean

solver-image:
	docker build \
	  --platform $(PLATFORM) \
	  --progress=plain \
	  --build-arg ELMER_DEBUG=$(ELMER_DEBUG) \
	  -f docker/Dockerfile \
	  -t $(SOLVER_IMAGE):latest \
	  .

clean:
	docker ps -a --filter "ancestor=$(SOLVER_IMAGE)" -q | xargs -r docker rm -f
	docker images $(SOLVER_IMAGE) -q | xargs -r docker rmi -f
