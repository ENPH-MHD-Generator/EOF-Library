# System packages and toolchain on top of eoflibrary/eof_elmer84_of6 (OpenFOAM 6
# and the Elmer source tree).
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

.PHONY: system

system:
	apt-get update
	# Debug tools, and METIS/ParMETIS (C libraries) for MUMPS
	apt-get install -y --no-install-recommends \
	  nano strace sudo wget ca-certificates software-properties-common \
	  libparmetis-dev libmetis-dev
	add-apt-repository -y ppa:ubuntu-toolchain-r/test
	apt-get update
	apt-get install -y --no-install-recommends gcc-9 g++-9 gfortran-9
	rm -rf /var/lib/apt/lists/*
	# Passwordless sudo for the openfoam user in interactive shells
	echo "openfoam ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/openfoam
	chmod 0440 /etc/sudoers.d/openfoam
	# Replace the EOF-Library copy that ships in the base image
	rm -rf $(EOF_HOME)
	mkdir -p $(EOF_HOME) $(EOF_LIB)
	sed -i '\|EOF-Library/etc/bashrc|d' /home/openfoam/.bashrc
	echo '. $(OPENFOAM_BASHRC)' >> /home/openfoam/.bashrc
