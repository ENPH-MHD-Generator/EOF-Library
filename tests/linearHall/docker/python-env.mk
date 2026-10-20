# Case preparation's Python environment, built inside the image's python stage:
# a uv-managed Python 3.11 with the locked dependencies of pyproject.toml
# (plasma_collisions from GitHub), and the BOLSIG+ distribution whose Phelps
# argon cross sections plasma_collisions tabulates. BOLSIG+'s terms of use are
# its authors'.

SHELL := /bin/bash
.SHELLFLAGS := -euo pipefail -c
.ONESHELL:

ENV_DIR := /opt/linear-hall/env

export UV_PYTHON_INSTALL_DIR := $(ENV_DIR)/python
export UV_PROJECT_ENVIRONMENT := $(ENV_DIR)/venv
export UV_LINK_MODE := copy
export UV_COMPILE_BYTECODE := 1
export PLASMA_COLLISIONS_BOLSIG_DIR := $(ENV_DIR)/bolsig

.PHONY: python-env

python-env:
	git config --global url."https://x-access-token:$$(cat /run/secrets/github_token)@github.com/".insteadOf "https://github.com/"
	uv python install 3.11
	uv sync --frozen --no-dev --no-cache
	$(UV_PROJECT_ENVIRONMENT)/bin/python -m plasma_collisions fetch
	$(UV_PROJECT_ENVIRONMENT)/bin/python -m plasma_collisions where
