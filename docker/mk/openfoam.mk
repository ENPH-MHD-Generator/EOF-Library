# OpenFOAM side: the EOF coupler library, then the OpenFOAM solvers that link it.
# Both install beneath $FOAM_USER_LIBBIN / $FOAM_USER_APPBIN (/home/openfoam/platforms).
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

# OpenFOAM's bashrc references unset variables and returns nonzero statuses
OPENFOAM_ENV := set +eu; source $(OPENFOAM_BASHRC); set -eu; export EOF_SRC=$(EOF_SRC)

.PHONY: coupler solvers

coupler:
	$(OPENFOAM_ENV)
	wclean $(EOF_SRC)/coupleElmer
	wmake -j $(JOBS) $(EOF_SRC)/coupleElmer

solvers:
	$(OPENFOAM_ENV)
	wclean $(EOF_HOME)/solvers/mdhLinearHall
	wmake -j $(JOBS) $(EOF_HOME)/solvers/mdhLinearHall
	test -x "$$FOAM_USER_APPBIN/mdhLinearHall"
