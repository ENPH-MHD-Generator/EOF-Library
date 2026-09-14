# Host-side build of the EOF-Library MHD images. The recipes that run inside the
# images live in docker/mk (solver image) and tests/linearHall (application).
#
#   make solver-image   EOF-Library solvers: OpenFOAM + Elmer + MHD modules
#   make image          solver image, then the linear Hall application on top
#   make clean          remove both images and their containers

SOLVER_IMAGE ?= eof-mhd-solvers
IMAGE_NAME   ?= mhd-sim
ELMER_DEBUG  ?= 0
PLATFORM     := linux/amd64

.PHONY: image solver-image clean

image: solver-image
	$(MAKE) -C tests/linearHall image IMAGE_NAME=$(IMAGE_NAME) SOLVER_IMAGE=$(SOLVER_IMAGE):latest

solver-image:
	docker build \
	  --platform $(PLATFORM) \
	  --progress=plain \
	  --build-arg ELMER_DEBUG=$(ELMER_DEBUG) \
	  -f docker/Dockerfile \
	  -t $(SOLVER_IMAGE):latest \
	  .

clean:
	$(MAKE) -C tests/linearHall clean IMAGE_NAME=$(IMAGE_NAME)
	docker ps -a --filter "ancestor=$(SOLVER_IMAGE)" -q | xargs -r docker rm -f
	docker images $(SOLVER_IMAGE) -q | xargs -r docker rmi -f
