# EOF-Library MHD solver image: OpenFOAM 6 with the EOF coupler and MHD solver,
# Elmer 8.4 (gcc 9, MUMPS, Hypre) with the EOF modules and MHDSolve.
#
#   base ─┬─ openfoam: coupler library → OpenFOAM solvers ─────┐
#         └─ elmer:    ScaLAPACK/MUMPS/Hypre → Elmer → extensions ┴─ solvers (final)
#
# The two branches share no build inputs, so editing one never rebuilds the
# other. Recipes live in docker/mk; each stage copies only the files it needs.
# Build with `make solver-image` from the repository root.

FROM eoflibrary/eof_elmer84_of6:latest AS base
USER root
COPY docker/mk/common.mk docker/mk/base.mk /build/mk/
RUN make -f /build/mk/base.mk system
WORKDIR /home/openfoam/EOF-Library


FROM base AS openfoam
COPY docker/mk/openfoam.mk /build/mk/
COPY libs/coupleElmer/ libs/coupleElmer/
COPY libs/commSplit/ libs/commSplit/
RUN make -f /build/mk/openfoam.mk coupler
COPY solvers/mdhLinearHall/ solvers/mdhLinearHall/
RUN make -f /build/mk/openfoam.mk solvers


FROM base AS elmer
COPY docker/mk/deps.mk /build/mk/
RUN make -f /build/mk/deps.mk deps
ARG ELMER_DEBUG=0
COPY docker/mk/elmer.mk /build/mk/
RUN --mount=type=cache,id=mhd-elmer-build-gcc9,target=/opt/elmerfem/build,sharing=locked \
    make -f /build/mk/elmer.mk elmer ELMER_DEBUG=${ELMER_DEBUG}
COPY docker/mk/extensions.mk /build/mk/
COPY libs/Elmer2OpenFOAM.F90 libs/OpenFOAM2Elmer.F90 libs/
COPY libs/solvers/MHDSolve/ libs/solvers/MHDSolve/
RUN make -f /build/mk/extensions.mk elmer-extensions


# Merge the branches: start from the Elmer stage and copy in what the OpenFOAM
# branch built. Both branches start FROM base, so OpenFOAM itself, MPI and the
# system libraries are already identical here; the OpenFOAM branch's only new
# outputs are under /home/openfoam/platforms ($FOAM_USER_LIBBIN: libElmer.so,
# $FOAM_USER_APPBIN: mdhLinearHall). Anything openfoam.mk installs elsewhere
# would be lost, so add a COPY for it (or build it in base). The merge starts
# from Elmer because its outputs are spread over /usr/local, /opt/mumps,
# /opt/hypre and /opt/eof/lib, while OpenFOAM's sit in one directory.
FROM elmer AS solvers
COPY --from=openfoam /home/openfoam/platforms/ /home/openfoam/platforms/
COPY --from=openfoam /home/openfoam/EOF-Library/solvers/ /home/openfoam/EOF-Library/solvers/

# Elmer finds the extension modules (Procedure = "MHDSolve" ...) on this path.
# One BLAS thread per MPI rank: the ranks already occupy the cores, and
# OpenBLAS's default of one thread per core in every rank made MUMPS 10-30x slower.
ENV LD_LIBRARY_PATH=/opt/eof/lib \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1
