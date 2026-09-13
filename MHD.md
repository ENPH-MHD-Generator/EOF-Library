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

All solver, compiler, Python, Elmer, OpenFOAM, Gmsh, and coupling dependencies
are installed in the image. The host does not need them.

## Build and enter the container

From the repository root, build the image once:

```sh
./mhd build
```

Open a container shell without rebuilding:

```sh
./mhd shell
```

`./mhd start` is a convenience command that builds and then opens the shell.
Pass `--debug` to `build` or `start` to compile Elmer's MHD module with runtime
checks and debug symbols:

```sh
./mhd start --debug
```

The image is named `mhd-sim:latest`. Docker build caching makes unchanged
rebuilds fast.

### Rebuild after changing solver source

Solver development uses the same reproducible image build instead of copying
individual files into a running container. After editing either the Elmer
Fortran module or the OpenFOAM C++ solver on the host, leave the current
container and run:

```sh
./mhd build
./mhd shell
```

The native builds are independent Docker stages:

- Changes beneath `libs/solvers/MHDSolve/` rebuild the Elmer MHD module but
  reuse the cached OpenFOAM solver and coupler.
- Changes beneath `solvers/mdhLinearHall/` rebuild the OpenFOAM solver but
  reuse the cached Elmer module and coupler.
- Changes to `libs/coupleElmer/`, `libs/commSplit/`,
  `libs/Elmer2OpenFOAM.F90`, or `libs/OpenFOAM2Elmer.F90` rebuild the shared
  coupler and both dependent solver stages.

The Elmer stage also keeps its configured MPI compiler tree in a BuildKit cache,
so changed Fortran sources can reuse previously compiled dependencies.

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

`mesh` is optional:

- `size_min` and `size_max` set optional global Gmsh edge-length bounds.
- `size_factor` scales Gmsh's characteristic lengths. Values below one refine
  the mesh and values above one coarsen it. The default is `1.0`.

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

`plasma` is optional and describes the alkali-seeded carrier gas. The defaults
are 1% potassium in argon:

```yaml
plasma:
  seed_mole_fraction: 0.01        # seed atoms per heavy particle, in (0, 1)
  seed_ionization_energy: 4.3407  # eV
  seed_gi_over_gn: 0.5            # ion / neutral statistical weight ratio
  seed_cross_section: 4.0e-18     # m^2, electron-seed momentum transfer
  carrier_cross_section: 1.0e-19  # m^2, electron-carrier momentum transfer
  reference_pressure: 101325      # Pa, absolute pressure where OpenFOAM p = 0
  sigma_min: 1.0e-2               # S/m
  sigma_max: 1.0e6                # S/m
  two_temperature: true           # Te from Joule heating vs. collisional loss
  carrier_molar_mass: 39.948      # g/mol
  seed_molar_mass: 39.098         # g/mol
  energy_loss_factor: 1.0         # delta; 1 = elastic losses only
  electron_temperature_max: 20000 # K
  electron_temperature_relaxation: 0.5
```

Only the seed ionizes. At every node Elmer solves the Saha equation for the
seed at the electron temperature `Te`, with heavy-particle densities from the
gas temperature, then computes the conductivity from electron-neutral
collisions with the carrier gas and the neutral seed. The Hall term uses the
resulting electron density, `1/(n_e e)`. Conductivity is clamped to
`[sigma_min, sigma_max]`, and the Hall parameter stays physical at clamped
nodes.

With `two_temperature: true`, `Te` comes from the electron energy balance
(Kerrebrock): Joule heating of the electrons, `J^2/sigma`, equals their
elastic collisional loss to heavy particles,
`3 delta n_e m_e k_B (Te - Tg) sum_s nu_s / M_s`. The heating uses the field
the electrons see, `E' = -grad(phi) + U x B`, with the Hall effect included:
`J^2/sigma = sigma (E'_par^2 + E'_perp^2 / (1 + beta^2))`, where
`beta = mu_e |B|`. Because conductivity depends
on `Te` and the current depends on conductivity, `Te` is updated every
nonlinear iteration of the current solver (under-relaxed by
`electron_temperature_relaxation`), and the solver only stops once `Te` has
also converged. `energy_loss_factor` scales the losses for inelastic or
radiative processes. With `two_temperature: false`, `Te` equals the gas
temperature.

The run writes `ionizationFraction` (`n_e / n_heavy`), `Te`, and
`elcond_elmer` as OpenFOAM fields, including initial values at `t = 0`. The
channel's initial temperature is the inlet temperature.

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
loses, `P_mech`, which should match `P_emf`. The flow is incompressible, so
extracted power appears as a pressure drop rather than a fall in gas
enthalpy.

`coupling` is optional and sets how often OpenFOAM re-solves the electrical
problem in Elmer:

```yaml
coupling:
  velocity_tolerance: 0.05        # max |U - U_sent| / max |U_sent|
  temperature_tolerance: 0.005    # max |T - T_sent| / T_sent
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

Use `--ranks 4` from about `size_factor` 0.2 down; on coarser meshes the fixed
start-up and coupling costs dominate and extra ranks do not help. The mesh is
Netgen-optimized after generation because sliver tetrahedra otherwise set the
time step for the whole mesh. Each Elmer update prints the electron
temperature convergence; nodes whose updates oscillate are damped
automatically, and the count appears as `damped`.

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
