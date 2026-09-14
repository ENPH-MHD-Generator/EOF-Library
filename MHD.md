# MHD solver guide

This guide is the complete user-facing workflow for the MHD solver in this
repository. Users configure one simulation in YAML; the tooling creates the
Elmer and OpenFOAM inputs, meshes, partitions, and run directory. Editing
`case.sif` or OpenFOAM dictionaries is not part of the supported workflow.

The pipeline deliberately supports one YAML file, one prepared case, and one
solver run. Parameter studies and multi-experiment manifests are not supported.
Use an external script or workflow system to invoke these commands repeatedly
when a sweep is needed.

## Prerequisites

- Docker with support for Linux AMD64 containers
- A checkout of this repository
- Read access to the private
  [`plasma_collisions`](https://github.com/ENPH-MHD-Generator/plasma_collisions)
  repository, through `gh auth login` or a `GITHUB_TOKEN` environment variable

All solver, compiler, Python, Elmer, OpenFOAM, Gmsh, and coupling dependencies
are installed in the image. The host does not need them.

## Build and enter the container

From the repository root, build the images once:

```sh
./mhd build
```

Open a container shell without rebuilding:

```sh
./mhd shell
```

`./mhd start` is a convenience command that builds and then opens the shell.
Pass `--debug` to `build` or `start` to compile Elmer with runtime checks and
debug symbols:

```sh
./mhd start --debug
```

Docker build caching makes unchanged rebuilds fast.

### How the images are built

There are two images, and the host only runs `docker build` (root `Makefile`,
`tests/linearHall/Makefile`); every build step inside an image is a `make`
target in a small recipe file:

- **`eof-mhd-solvers:latest`**, the EOF-Library solver image
  (`docker/Dockerfile`, recipes in `docker/mk/`). On top of
  `eoflibrary/eof_elmer84_of6` (OpenFOAM 6 and the Elmer source) it builds two
  independent branches and merges them:
  - OpenFOAM: the EOF coupler library (`libs/coupleElmer`), then the OpenFOAM
    solvers (`solvers/mdhLinearHall`) that link it (`openfoam.mk`).
  - Elmer: ScaLAPACK, MUMPS and Hypre (`deps.mk`), Elmer itself (`elmer.mk`),
    then the Elmer extensions, the EOF coupler modules and MHDSolve
    (`extensions.mk`), installed to `/opt/eof/lib`.
- **`mhd-sim:latest`**, the linear Hall application (`tests/linearHall/Dockerfile`,
  built from that folder alone): case preparation's Python environment
  (`tests/linearHall/docker/python-env.mk`, with `plasma_collisions` and
  BOLSIG+), the case template, and the `mhd` command, `FROM` the solver image.

`./mhd` at the repository root forwards to `tests/linearHall/mhd`, which builds
the solver image from this checkout and then the application. The application
folder has no other dependency on the rest of the repository: moved elsewhere,
its `mhd` builds against an existing `eof-mhd-solvers` image (or one named by
`SOLVER_IMAGE`), with `MHD_WORKSPACE` choosing where `experiments/` and `out/`
live. `make solver-image` and `make image` in the repository root build the
images without the launcher.

### Rebuild after changing solver source

Solver development uses the same reproducible image build instead of copying
individual files into a running container. After editing on the host, leave the
current container and run:

```sh
./mhd build
./mhd shell
```

What rebuilds:

- `libs/solvers/MHDSolve/`, `libs/Elmer2OpenFOAM.F90`, `libs/OpenFOAM2Elmer.F90`:
  the Elmer extensions only.
- `solvers/mdhLinearHall/`: the OpenFOAM solver only.
- `libs/coupleElmer/`, `libs/commSplit/`: the coupler and the OpenFOAM solver.
- `docker/mk/*.mk`: the stage using that recipe and those after it in its
  branch (for example `elmer.mk` rebuilds Elmer, incrementally, and the
  extensions, but not the third-party libraries or the OpenFOAM branch).
- Anything in `tests/linearHall/`: the application image only;
  `pyproject.toml` or `uv.lock` also rebuild its Python environment. To move to
  a newer `plasma_collisions`, run `uv lock --upgrade-package plasma-collisions`
  in `tests/linearHall`.

Elmer and its Fortran dependencies are built with gcc/gfortran 9 (Ubuntu
toolchain PPA); OpenFOAM, the OpenFOAM coupler and Open MPI 1.10 keep the system
gcc 5. Elmer only uses `mpif.h`, so Open MPI's wrappers are pointed at gfortran 9
with `OMPI_FC` for those builds. ScaLAPACK 2.1.0, MUMPS 5.6.2 (parallel direct
solver) and Hypre 2.15.1 are built from checksummed source tarballs with the
same compiler: Fortran libraries that share derived types with Elmer must use
one compiler, since gfortran 8 changed the array descriptor ABI. Hypre is
available for future symmetric problems; BoomerAMG does not converge on the
non-symmetric, penalty-coupled potential equation. BLAS runs one thread per MPI
rank (`OPENBLAS_NUM_THREADS=1`).

The Elmer extensions are always compiled after Elmer, against the installed
version: Elmer's data types depend on its build options, so modules built
against a different Elmer would corrupt memory. Elmer's build directory is a
BuildKit cache, so changing its options (including `--debug`) recompiles Elmer
incrementally.

Prepared cases and results remain in host `out/`, so replacing the disposable
container does not remove them. Once back inside, rerun a prepared case with
`mhd run NAME`, or rerun `mhd prepare ... --force` first if configuration or
case-generation code changed. Use `./mhd build --debug` when developing the
Elmer module with runtime checks.

## Configure one case

Put experiment files in the host repository's `experiments/` directory. The
container mounts that directory read-only at `/experiments`. Start by copying
[`experiments/linear-hall.yaml`](experiments/linear-hall.yaml).

```yaml
schema_version: 1

channel:
  length: 0.200
  height: 0.050
  width: 0.050
  wall_thickness: 0.005

mesh:
  size_min: null
  size_max: null
  size_factor: 0.2

electrodes:
  length: 0.010
  pairs:
    - {x_center: 0.040, resistance: 1.0}
    - {x_center: 0.080, resistance: 1.0}
    - {x_center: 0.120, resistance: 1.0}
    - {x_center: 0.160, resistance: 1.0}

physics:
  B_field: [0, 0, 1.5]
  inlet_velocity: [391, 0, 0]
  inlet_temperature: 1500
```

All geometry values are in metres. `B_field` is in tesla,
`inlet_velocity` is in metres per second, `inlet_temperature` is in kelvin, and
electrode resistance is in ohms.

### Schema reference

`schema_version` must be `1`.

`channel` defines the rectangular flow channel and is required:

- `length`, `height`, and `width` must be positive.
- `wall_thickness` is the positive thickness used for the electrode and
  insulating shell geometry.

`mesh` is optional. `type` selects the mesh:

- `tetrahedral` (default): unstructured tetrahedra.
  - `size_min` and `size_max` set optional global Gmsh edge-length bounds.
  - `size_factor` scales Gmsh's characteristic lengths. Values below one refine
    the mesh and values above one coarsen it. The default is `1.0`.
- `structured`: graded hexahedra that resolve the cold thermal boundary layer
  on the electrode walls, refined streamwise at the electrode edges.
  - `electrode_wall_cell_size` (default 1.5e-5 m) is the first cell on the walls
    that carry the electrodes (y = 0 and y = height).
  - `side_wall_cell_size` (default: `wall_cell_size`) is the first cell on the
    insulating side walls (z = 0 and z = width).
  - `wall_cell_size` (default 5e-4 m) is the first cell on walls without a
    specific size.
  - `cell_size` (default 0.003 m) is the core cell size across the channel, and
    `streamwise_cell_size` (default 0.005 m) the core size along it.
  - `electrode_edge_cell_size` (default 5e-4 m) is the streamwise size at the
    electrode edges.
  - `growth_rate` (default 1.2) is the largest size ratio of neighbouring cells.

  Keys of the other mesh type are rejected. Setting `electrode_wall_cell_size`
  or `side_wall_cell_size` to `null` uses `wall_cell_size`, and
  `streamwise_cell_size: null` uses `cell_size`.

  Structured meshes need the electron energy equation (the default) or
  `two_temperature: false`: with the local electron energy balance, their wall
  cells are far smaller than the ~1 mm electron energy relaxation length and
  the balance runs away at the electrode edges. The grid is a tensor product,
  so a fine wall layer runs the whole channel length and a fine electrode-edge
  slice spans the whole cross-section. The defaults give about 70k nodes for an
  80 mm channel with one electrode pair and about 230k for the 200 mm channel
  with four; see "Resolving the cold wall layer" below for their accuracy and
  cost.

`electrodes` is required. Every electrode pair must have the same positive
`length`; pair-specific lengths are rejected. Choose exactly one placement
form.

Explicit placement supports a resistance for each pair:

```yaml
electrodes:
  length: 0.010
  pairs:
    - {x_center: 0.050, resistance: 1.0}
    - {x_center: 0.150, resistance: 2.0}
```

Automatic placement distributes a requested count evenly along the channel and
uses one resistance for all pairs:

```yaml
electrodes:
  length: 0.010
  count: 4
  resistance: 1.0
```

Pairs may not overlap or extend beyond the channel. Resistance must be
non-negative.

`physics` is optional:

- `B_field` is a three-component magnetic flux-density vector.
- `inlet_velocity` is a three-component inlet velocity vector.
- `inlet_temperature` must be positive.
- `insulator_wall_temperature` and `electrode_wall_temperature` (default 300 K)
  fix the temperature of the insulating walls and the electrodes, so heat is
  lost to them. `null` makes that surface adiabatic. The defaults assume
  uncooled walls starting at room temperature: over a run of under 10 s a
  copper electrode surface warms about 10 K and a ceramic one tens to a few
  hundred K. Elmer sees the cold wall through the temperature interpolation,
  so the conductivity drops in the cold layer next to the walls.
- `outlet_pressure` (default 101325 Pa) is the absolute static pressure at the
  outlet.

The gas is argon, solved as a compressible, laminar, calorically perfect ideal
gas (`rhoPimpleFoam`; `constant/thermophysicalProperties`): `Cp` = 520.3
J/(kg K), Sutherland viscosity fitted to argon, and conductivity from the
Eucken relation (Prandtl number 2/3). Density therefore follows the
temperature, which matters in the cold wall layer (about 8x denser than a
2500 K core), and the extracted electrical power leaves the gas as enthalpy.
The initial velocity is potential flow, and the initial pressure is uniform at
`outlet_pressure`. Starting from that state launches pressure waves while the
boundary layers form; on the tetrahedral meshes the cell pressure then stays
within about 10% of the mean, with the extremes in single corner cells. The
log prints the density, pressure, temperature and Mach number range each step.
The schemes in `system/fvSchemes` are chosen for robustness on the
tetrahedral meshes: Euler time stepping and bounded (`limitedLinear`)
convection of enthalpy and kinetic energy. With
`backward` and `linearUpwind`/`linear` energy convection, corner cells heated
spuriously by ~1000 K and the run diverged within 4e-5 s even without a
magnetic field.

`plasma` is optional and describes argon seeded with potassium:

```yaml
plasma:
  seed_mole_fraction: 0.01          # potassium atoms per heavy particle, in (0, 1)
  electron_transport_model: drifting  # or lorentz
  potassium_elastic_scale: 1.0      # multiplier on the e-K momentum cross section
  excited_state_temperature: electron  # or gas
  ion_reduced_mobility: 2.43e-4     # m^2/(V s) at 2.6868e25 m^-3 (K+ in Ar)
  sigma_min: 1.0e-2                 # S/m
  sigma_max: 1.0e6                  # S/m
  two_temperature: true             # Te from Joule heating vs. collisional loss
  electron_energy_transport: true   # electron energy equation (false: local balance)
  electron_temperature_max: 20000   # K
  electron_temperature_relaxation: 0.5
```

Only the seed ionizes: at every node Elmer solves the Saha equation for
potassium at the electron temperature `Te`, with heavy-particle densities from
the gas temperature and pressure.

Electron collision data come from the
[`plasma_collisions`](https://github.com/ENPH-MHD-Generator/plasma_collisions)
package, for a Maxwellian electron energy distribution at `Te`: argon's Phelps
cross sections (including its Ramsauer minimum) and potassium cross sections
from the literature. `mhd prepare` tabulates what depends on `Te` alone into
`electron_collisions.dat` in the case (provenance in
`electron_collisions.json`), and MHDSolve reads it once and combines it with
the local state:

- **Momentum transfer and Coulomb collisions.** Electron-neutral collisions
  come from the tables; electron-ion (Coulomb) collisions use the NRL
  collision frequency and Coulomb logarithm at the local `n_e` and `Te`. In the
  channel core they are about half of all electron collisions, and most of them
  near the electrodes.
- **`electron_transport_model`.** A Maxwellian does not fix how electrons of
  different speeds share the drift, and argon's Ramsauer minimum makes that
  matter. `drifting` (the default) assumes strong electron-electron collisions:
  collision frequencies add, giving one collision frequency and the familiar
  tensor with `beta = mu_e |B|`; it is a lower bound on the conductivity.
  `lorentz` assumes none: collision frequencies add inside the velocity
  average, the Pedersen and Hall mobilities are tabulated against `Te`,
  `n_e ln(Lambda)/N` and `omega_ce/N` (a 9 MB table, within 3% of the
  library), and it is an upper bound. Where electron-electron collisions are
  comparable to the others, as here, the truth lies between them. At the
  benchmark conditions they differ by ~15% in the Pedersen conductivity and ~2%
  in load power.
- **Ion slip.** Ions carry current too, with mobility
  `ion_reduced_mobility x N0/N` (polarization collisions; the default is the
  Langevin value for K+ in argon). Their Hall current runs opposite to the
  electrons', which reduces the effective Hall parameter when `beta_e beta_i`
  approaches 1: negligible at 1 atm (~1%), significant below ~0.1 atm.
- **Conductivity tensor.** The current uses the full tensor,
  `J = sigma_0 (b.E') b + sigma_P E'_perp + sigma_H b x E'` with `b = B/|B|`,
  from electrons and ions. For a single collision frequency this is exactly the
  generalized Ohm's law `E' = J/sigma + J x B/(n_e e)`. The parallel
  conductivity is clamped to `[sigma_min, sigma_max]`, scaling the Pedersen and
  Hall components with it so the Hall parameter stays physical.
- **Energy losses.** Elastic recoil to neutrals and ions, from the tables.
  `excited_state_temperature` sets the K(4p) population that returns
  excitation energy through superelastic collisions: `electron` (default)
  assumes electron collisions keep it at `Te`, so excitation and de-excitation
  balance and there is no net inelastic loss (resonance radiation is trapped in
  the channel); `gas` assumes the excited atoms are quenched to the gas
  temperature, which makes resonance excitation the dominant electron energy
  loss, about 1000x the elastic loss at `Te` = 3500 K.
- **`electron_wall_heat_transmission`** (default 0: walls adiabatic for the
  electrons). The electron energy carried into walls and electrodes through a
  sheath: electrons arriving at the Bohm flux `0.61 n_e sqrt(k_B Te / M_K+)`
  each deposit this many `k_B Te`. About 6.7 for a floating wall with K+ ions
  (`2 + ln(M/(2 pi m_e))/2`). Current-carrying electrodes are treated the same
  way, which is an approximation.
- **`potassium_elastic_scale`.** The e-K momentum-transfer cross section is
  the least certain input (+-30%, from caesium data); scale it for
  sensitivity studies.

Building the application image fetches `plasma_collisions` from its private
GitHub repository (see `tests/linearHall/pyproject.toml` and `uv.lock`), using a
GitHub token passed to Docker as a build secret: the GitHub CLI login
(`gh auth login`) or `GITHUB_TOKEN`. The BOLSIG+ distribution, whose Phelps
argon cross sections the tables use, is downloaded into the image during the
build; its terms of use are its authors', and publications must cite Hagelaar
and Pitchford, Plasma Sources Sci. Technol. 14, 722 (2005).

With `two_temperature: true`, `Te` rises above the gas temperature where Joule
heating of the electrons, `J_e.E'`, exceeds their collisional energy losses.
With `two_temperature: false`, `Te` equals the gas temperature.

With `electron_energy_transport: true` (the default), Elmer solves the steady
electron energy equation (Solver 11 in `case.sif`):

```text
(5/2) k_B Gamma_e . grad(Te) + ((5/2) k_B Te + chi) div(n_e U) - div(K_e grad(Te))
    = J_e . E' - L_elastic - L_inelastic
```

- `Gamma_e = n_e U - J_e/e` is the electron flux: electron enthalpy is carried
  by the gas flow and by the electron current (from cathode to anode). `J_e` is
  the electron part of the current, with the electron conductivity tensor.
- `K_e = (5/2) (k_B^2 Te / e^2) sigma_e` is the electron heat conduction
  tensor, anisotropic like the conductivity: along `B` it spreads the heating
  over the relaxation length (about 1 mm at 1 atm, growing as 1/p), and across
  `B` it is weaker by about `1 + beta^2`. It keeps `Te` and the conductivity
  finite at the electrode edges and lets hot electrons reach into the cold wall
  layer.
- `chi` is the seed ionization energy. With Saha equilibrium `n_e` follows
  `Te`, so ionizing the seed along the flow takes `chi` per electron from the
  electrons (recombination returns it). This makes the electron temperature
  respond to the heating over a few centimetres downstream. Finite-rate ionization
  (for example from an inductively coupled plasma source) would replace the
  Saha relation and this term with an electron continuity equation.
- The steady form is valid because electron energy relaxes in ~1e-7 s, far
  faster than the flow changes.
- Incoming gas is in equilibrium, `Te = Tg`, at the inlet. The equation is
  dominated by convection there, and without that inflow value `Te` drifted
  below the gas temperature and flickered over the first few centimetres.
  Other boundaries have no conductive electron heat flux; walls and electrodes
  cool the electrons through the gas temperature they relax to.
- The heating and the electron drift use the physical conductivity
  `e n_e mu_e`, not the value clamped to `sigma_min`, so cold gas does not
  heat its vanishing electron population.

The equation is linearized in `n_e(Te)` (one Newton step per nonlinear
iteration), stabilized with SUPG, and solved inside the current solver's
nonlinear iteration because `Te`, the conductivity and the current are tightly
coupled. Each iteration's `Te` update is under-relaxed by
`electron_temperature_relaxation`, and the current solver stops once both the
potential and `Te` have converged. The Elmer log prints the electron energy
budget, `Electron energy: Joule ... elastic loss ... inelastic loss ...
ionization ...`, whose terms should balance, and each update prints the largest
electron Hall parameter and share of electron-ion collisions.

With `electron_energy_transport: false`, `Te` comes from the local balance of
heating and the same losses at each node, which ignores all transport. It is
cheaper but runs away at electrode edges on fine meshes.

The run writes `ionizationFraction` (`n_e / n_heavy`), `Te`, `elcond_elmer`
(the parallel conductivity), `pedersenConductivity` and `hallParameter` (the
electron Hall parameter, Hall over Pedersen mobility) as OpenFOAM fields,
including initial values at `t = 0`. The channel's initial temperature is the
inlet temperature. The Lorentz damping in the momentum equation uses the
Pedersen conductivity.

Elmer receives the velocity as cell values (`interpolationSchemes` in
`system/fvSchemes` uses `cell` for `Ux/Uy/Uz`, not `cellPoint`). `cellPoint`
blends in vertex values, which are exactly zero on no-slip walls, so every
Elmer node on a wall would lose its `U x B` EMF across a whole cell. Because
the boundary layer is far thinner than a cell here, that removed most of the
generated power. The Elmer log prints the volume-averaged velocity it sees:
it should be close to the inlet velocity, and a much lower value means the
EMF is being damped this way.

The current and Joule heating are transferred per element rather than
interpolated from nodes, so the peaks at the electrode edges reach OpenFOAM
undiminished.

The full Lorentz force `J x B` acts on the bulk gas. Electrons and ions pass
their momentum to the neutrals by collisions within nanoseconds, and ion slip
(`beta_i ~ 1e-3`) is negligible. The force uses the current from the previous
coupling step; a semi-implicit drag, `-sigma |B|^2 (U - U_sent)`, stabilizes
that lag and vanishes at convergence. The Elmer log prints a power balance,
`P_emf = int J.(U x B) dV` against Joule heating plus the power delivered to
the electrode loads, and the OpenFOAM log prints the mechanical power the flow
loses, `P_mech`, which should match `P_emf`. The energy equation receives the
electromagnetic power `J.E = J^2/sigma + U.(J x B)`: Joule heating less the
extracted mechanical power. It uses the same force as the momentum equation,
so the gas loses exactly the power delivered to the loads as enthalpy.

`coupling` is optional and sets how often OpenFOAM re-solves the electrical
problem in Elmer:

```yaml
coupling:
  velocity_tolerance: 0.05        # max |U - U_sent| / max |U_sent|
  temperature_tolerance: 0.005    # Joule-power-weighted RMS |T - T_sent| / T_sent
  pressure_tolerance: 0.05        # max |p - p_sent| / absolute pressure
  max_steps_between_updates: 250  # 0 = no limit
```

The electrical problem is quasi-static: the current depends only on the
instantaneous velocity, temperature and pressure, so skipping an Elmer update
only means using slightly stale inputs. Elmer is re-solved once any input has
changed by more than its relative tolerance since the last update, after
`max_steps_between_updates` steps, and always on the final step. All zero (the
default) updates every time step. Between updates the Lorentz force is held,
with the implicit drag term correcting it for the velocity change.

Conductivity is exponential in temperature (about 1% per 2.5 K at 2500 K), so
keep `temperature_tolerance` roughly ten times tighter than the velocity
tolerance. To check a setting, read the OpenFOAM log: each update prints the
input changes that triggered it and `change since last update`, the relative
jump in `J x B` when it is refreshed, which is the error the skipped updates
carried. Keep it to a few percent, and compare integral results (`P_emf`,
electrode currents) against a run with all tolerances at zero.

`numerics` is optional:

```yaml
numerics:
  linear_solver: auto   # auto, iterative, or mumps
```

`linear_solver` selects how Elmer solves the potential equation: `iterative`
(ILU-preconditioned GCR) or `mumps` (parallel sparse direct). `auto` uses
iterative on tetrahedral meshes, where it was about 20% faster at
`size_factor` 0.15, and MUMPS on structured meshes. Each iterative solve starts
from the previous solution, so once a run is going GCR needs only ~50-150
iterations, while MUMPS refactors the matrix every time. On the 72k-node
structured test mesh `iterative` was therefore 12% faster overall than MUMPS
(potential solves 57 s against 81 s), but on the 230k-node full-length
structured mesh its first solves, from a zero potential on thin wall cells,
took ~215 s each against ~10 s for MUMPS. Structured meshes with `iterative`
use row scaling and a tolerance of 1e-5, which matched MUMPS to 0.2% in load
power (1e-4 left 1.6%).

Unknown keys and invalid values fail during validation instead of being
silently ignored.

## Prepare one case

List every `.yaml` or `.yml` experiment currently visible in the host-mounted
`experiments/` directory:

```sh
mhd prepare --list
```

Inside the container, one command validates the YAML and produces an
execution-ready case:

```sh
mhd prepare /experiments/linear-hall.yaml
```

The default case name is the YAML filename without `.yaml`, so this writes
`/runs/linear-hall`. The output root is fixed: neither the YAML nor the command
can redirect a case outside `/runs`.

Preparation performs the following work:

1. Validates and resolves the YAML configuration.
2. Renders `case.sif` and all generated OpenFOAM initial fields.
3. Copies the maintained base `constant/` and `system/` inputs.
4. Reuses or generates the Gmsh mesh.
5. Converts the mesh for OpenFOAM and Elmer.
6. Runs `potentialFoam`, partitions both solvers with the same MPI rank count,
   validates the result, and writes a preparation marker.
7. Atomically exposes the completed case at `/runs/<name>`; a failed
   preparation does not leave a partially prepared named case.

By default, preparation reuses the checked-in baseline mesh. This is safe when
channel dimensions, mesh settings, electrode count and centers, and the shared
electrode length match the baseline. Changes limited to field values or
electrode resistance do not require a new mesh.

Use `--mesh` after changing any mesh-affecting setting:

```sh
mhd prepare /experiments/new-geometry.yaml --mesh
```

Useful options are:

- `--name NAME` chooses the directory name beneath `/runs`.
- `--ranks N` partitions each solver into `N` MPI ranks; the default is `2`.
- `--force` atomically replaces an existing case with the same name.
- `--dry-run` validates the input and prints the planned files and commands
  without writing a case.

Case names may contain letters, numbers, `.`, `_`, and `-`, and must begin with
a letter or number.

## Run one case

List prepared cases that pass the same marker, rank, and mesh preflight checks
used by `mhd run`:

```sh
mhd run --list
```

Inside the same or a later container, execute the prepared case with one
command:

```sh
mhd run linear-hall
```

The command only accepts a prepared case name, not an arbitrary path. It reads
the rank count recorded during preparation, verifies the preparation marker and
required solver meshes, launches the coupled OpenFOAM/Elmer MPI process, and
then runs `reconstructPar` and `foamToVTK`.

OpenFOAM and Elmer take turns computing, so both solvers can share the same
cores: `mhd run` always allows oversubscription and has waiting MPI ranks yield
the CPU. `--cores C` confines the whole run to `C` cores, for example to match a
cluster allocation:

```sh
mhd run linear-hall --cores 4
```

With `--ranks 4`, running all 8 processes on 4 cores took about 20% longer than
spreading them over 10, for 60% fewer cores. After the simulation,
`reconstructPar` and `foamToVTK` run concurrently over the time directories,
one process per core (or per rank without `--cores`); the output is identical
to a serial run and each process logs to `log.<tool>.<n>` in the case.

To retain only the raw parallel output and skip reconstruction and VTK export:

```sh
mhd run linear-hall --no-postprocess
```

Validate and print the execution command without launching the solvers:

```sh
mhd run linear-hall --dry-run
```

### Run time and mesh size

Cost grows roughly as the fourth power of mesh refinement: halving the cell
size gives about 8x the cells and, through the Courant limit, about 2x the
time steps. Measured for 1e-4 s of simulated time:

| `size_factor` | Elements | 2 ranks | 4 ranks |
|---|---|---|---|
| 0.25 | ~18k | 20 s | |
| 0.2 | ~33k | 27 s | 20 s |
| 0.15 | ~77k | 100 s | 73 s |

These timings predate the fixed wall temperatures, the one-point tetrahedral
quadrature in Elmer and `maxCo 0.8`; the latter two cut run time by roughly
25% and 10%. They also predate the compressible solver and the electron energy
equation. For 6e-5 s at `size_factor` 0.15 with 4 ranks: incompressible flow
with the local balance took 37 s, compressible flow with the local balance
26 s, and compressible flow with the electron energy equation 22 s. The
compressible runs trigger fewer Elmer updates, and the electron energy
equation needs fewer nonlinear iterations per update than the local balance
(1.6 against 3.0); its own assembly and solve take about 1.5 s.

Near-electrode convergence was checked on an 80 mm channel with one electrode
pair at 4e-5 s (load power averaged over the last half). Tetrahedra of 4, 3 and
2 mm gave 134, 145 and 149 W; graded hexahedra (coarse, default and half the
default wall cell) gave 154, 153 and 152 W, with electrode currents within
0.6%. Integral results converge; the peak `Te` and conductivity at the
electrode edges keep rising slowly with refinement (the current concentration
there is singular), so read local maxima as mesh-dependent. `size_factor` is
relative to the domain size, so use `size_max` for absolute cell sizes.

### Resolving the cold wall layer

The gas next to the cooled walls is far colder than the core, and its low
conductivity is a classic source of large electrode voltage losses. The layer is
thin: its thickness grows with the transit time as `sqrt(alpha x / u)`, about
0.1 mm at an electrode 30 mm from the inlet after one transit and 0.5 mm after a
few. Tetrahedral meshes and 250 um wall cells cannot resolve it, and it only
forms once the gas has passed the electrode, so run for at least about twice
`x_electrode / u`.

With the structured defaults (15 um first cells on the electrode walls), an
80 mm channel with one electrode pair was followed to 1.5e-4 s:

- The gas temperature at the electrode centre rises from ~470 K in the first
  cell to the 2500 K core over ~0.5 mm, still slowly thickening.
- The electrons stay hot across it (Te ~2000 K in the first cell): Joule
  heating is intense where the conductivity is low, and with Saha ionization at
  Te the conductivity stays at ~10 S/m, far above `sigma_min`.
- The voltage drop across the layer (the wall potential against the core
  profile extrapolated to the wall) is 0.65-0.7 V at the cathode and
  0.25-0.3 V at the anode, against 12.4 V across the core. Halving the first
  cell to 7.5 um changed it by ~0.05 V; lowering `sigma_min` to 1e-6 S/m
  changed nothing.
- `electron_wall_heat_transmission: 6.7` (sheath energy loss) cools the
  electrons in the first ~30 um (Te 1720 K, conductivity 2 S/m at the wall) and
  takes 14% of the electron heating, but moves the voltage drop by only a few
  hundredths of a volt.

These drops are small because ionization is instantaneous (Saha at Te) and
charged particles are not lost to the walls. In the cold layer three-body
recombination takes ~1e-4 s while ambipolar diffusion reaches ~1 mm in that
time, so real walls should deplete the electrons across the whole layer. Finite-
rate ionization with ambipolar diffusion (and, at the cathode, a limit on
electron emission) is the physics that sets large boundary-layer losses; the
structured mesh resolves the layer for it.

Cost: that 80 mm case (72k nodes, 4 ranks) takes about 170 s per 5e-5 s of
simulated time, of which Elmer is ~75%, mostly the MUMPS potential solve
with MUMPS (~1.7 s per solve); `linear_solver: iterative` took about 150 s.
The electron energy system stays iterative: MUMPS was 4x slower for it. The full 200 mm channel (230k nodes) is several times that: MUMPS takes ~10 s
per potential solve there, so a resolved full-length run costs roughly 10-15 s
per time step on 4 ranks. It fits in memory (6 GB peak).

Use `--ranks 4` from about `size_factor` 0.2 down; on coarser meshes the fixed
start-up and coupling costs dominate and extra ranks do not help. The mesh is
Netgen-optimized after generation because sliver tetrahedra otherwise set the
time step for the whole mesh. Each Elmer update prints the electron
temperature convergence; nodes whose updates oscillate are damped
automatically, and the count appears as `damped`. Convergence uses the
volume-weighted RMS change of `Te`, so the few cells at the singular current
concentration at electrode edges no longer force every update to the
iteration limit on fine meshes; the largest nodal change is still logged.

## Where cases and results live

`./mhd shell` creates two host directories and bind-mounts them:

| Host path | Container path | Access | Purpose |
| --- | --- | --- | --- |
| `experiments/` | `/experiments` | read-only | User YAML inputs |
| `out/` | `/runs` | read/write | Prepared cases and solver results |

Consequently, `/runs/linear-hall` in the container is the same directory as
`out/linear-hall` on the host. Solver output appears there immediately and
survives when the disposable container exits. A separate copy or sync step is
not required.

These are live bind mounts. Adding or editing a host file beneath
`experiments/` is immediately reflected by `mhd prepare --list`, and a case
created beneath host `out/` is immediately reflected by `mhd run --list`. An
image rebuild is only required for source code or image dependency changes.

A prepared case includes the resolved YAML, a generation manifest, rendered
Elmer and OpenFOAM inputs, converted and partitioned meshes, and a hidden
`.mhd-prepared.json` execution marker. Running the solver adds its normal time
directories, logs/data, reconstructed fields, and `VTK/` output to the same
named directory.

`out/` is ignored by Git. Copy results elsewhere if they need long-term
archival outside the working tree.

## Complete example

From the host:

```sh
./mhd build
./mhd shell
```

Then inside the container:

```sh
mhd prepare /experiments/linear-hall.yaml
mhd run linear-hall
```

Inspect the result on the host at `out/linear-hall/`.

## Common failures

- **The case already exists:** choose another `--name`, or deliberately replace
  it with `--force`.
- **The configuration does not match the reusable mesh:** rerun preparation
  with `--mesh`.
- **A case is not prepared:** use `mhd prepare` first; `mhd run` will not run a
  hand-assembled or partial directory.
- **The image does not exist:** run `./mhd build` before `./mhd shell`.
- **Architecture warnings or failures:** ensure Docker can run `linux/amd64`
  images on the host.
