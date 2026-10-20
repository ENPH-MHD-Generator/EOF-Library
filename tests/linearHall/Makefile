# Host-side build of the linear Hall application image, on top of an existing
# EOF-Library solver image (build that with `make solver-image` in EOF-Library,
# or use ./mhd build).
#
#   make image   application image $(IMAGE_NAME):latest from $(SOLVER_IMAGE)
#   make clean   remove the application image and its containers
#
# plasma_collisions (a private GitHub repository) is fetched during the build
# with GITHUB_TOKEN if set, otherwise the GitHub CLI login.

IMAGE_NAME   ?= mhd-sim
SOLVER_IMAGE ?= eof-mhd-solvers:latest
PLATFORM     := linux/amd64

.PHONY: image clean

image: export GITHUB_TOKEN ?= $(shell gh auth token 2>/dev/null)
image:
	@test -n "$$GITHUB_TOKEN" || { echo "No GitHub token: run 'gh auth login' or set GITHUB_TOKEN (needs read access to ENPH-MHD-Generator/plasma_collisions)" >&2; exit 1; }
	@docker image inspect $(SOLVER_IMAGE) > /dev/null 2>&1 || { echo "Solver image $(SOLVER_IMAGE) not found: build it first (make solver-image in EOF-Library)" >&2; exit 1; }
	docker build \
	  --platform $(PLATFORM) \
	  --progress=plain \
	  --secret id=github_token,env=GITHUB_TOKEN \
	  --build-arg SOLVER_IMAGE=$(SOLVER_IMAGE) \
	  -t $(IMAGE_NAME):latest \
	  .

clean:
	docker ps -a --filter "ancestor=$(IMAGE_NAME)" -q | xargs -r docker rm -f
	docker images $(IMAGE_NAME) -q | xargs -r docker rmi -f
