# Elmer extensions: the Elmer side of the EOF coupler and the MHD solver module,
# installed to $(EOF_LIB). MHDSolve links the plasma library (plasma.mk). They must be compiled against the installed Elmer:
# its derived types change with build options (HAVE_MUMPS adds matrix fields),
# so modules built against another Elmer corrupt memory. Elmer finds
# `Procedure = "MHDSolve"` etc. through LD_LIBRARY_PATH.
include $(dir $(lastword $(MAKEFILE_LIST)))common.mk

MODULE_DIR := /tmp/eof-modules

.PHONY: elmer-extensions

elmer-extensions:
	export $(GCC9_ENV)
	mkdir -p $(EOF_LIB) $(MODULE_DIR)
	elmerf90 -o $(EOF_LIB)/Elmer2OpenFOAM.so -J $(MODULE_DIR) $(EOF_SRC)/Elmer2OpenFOAM.F90
	elmerf90 -o $(EOF_LIB)/OpenFOAM2Elmer.so -J $(MODULE_DIR) $(EOF_SRC)/OpenFOAM2Elmer.F90
	# MHDPlasma.F90 declares the plasma library (plasma/, built by plasma.mk),
	# whose static archive goes in after the objects that call it
	elmerf90 -o $(EOF_LIB)/MHDSolve.so -J $(MODULE_DIR) \
	  $(EOF_SRC)/solvers/MHDSolve/MHDUtils.F90 $(EOF_SRC)/solvers/MHDSolve/MHDPlasma.F90 \
	  $(EOF_SRC)/solvers/MHDSolve/MHDSolve.F90 \
	  $(PLASMA_LIB) -lpthread -ldl -lm
	rm -rf $(MODULE_DIR)
