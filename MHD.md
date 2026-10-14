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
schema_version: 2

channel:
  length: 0.200
  height: 0.050
  width: 0.050
  wall_thickness: 0.005

mesh:
  target_element_size: 0.005

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

openfoam:
  files:
    constant/transportProperties:
      entries:
        transportModel: Newtonian
        rho: 0.325
        nu: 2.13e-4
        cp: 520.0
        kappa: 0.0866
        elcond_melt: 3.0
        maxRelDiff: 0
    system/controlDict:
      entries:
        application: mdhLinearHall
        endTime: 1.0e-8
        adjustTimeStep: true
        maxCo: 0.002
```

All geometry values are in metres. `B_field` is in tesla,
`inlet_velocity` is in metres per second, `inlet_temperature` is in kelvin, and
electrode resistance is in ohms.

### Schema reference

New configurations use `schema_version: 2`. Version 1 remains loadable for
existing configurations, but does not support data-driven OpenFOAM files.

The versioned Pydantic models in
`tests/linearHall/caseprep/models.py` are the authoritative schema. They reject
unknown keys, validate types and ranges, and enforce physical constraints that
span configuration sections. The loader selects the matching model using
`schema_version`, which allows future versions to coexist without silently
changing the meaning of existing files.

The generated JSON Schema at `schemas/mhd-experiment-v2.schema.json` provides
portable editor validation and completion. The example YAML declares this
schema with a `yaml-language-server` comment. After changing a Pydantic model,
regenerate and verify the artifact from a Python environment containing
`requirements-caseprep.txt`:

```sh
make schema
make check-schema
```

Both versioned schemas are generated and checked. For local development, use
a virtual environment and run the config and pipeline regression suite:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-caseprep.txt
make test-caseprep PYTHON=.venv/bin/python
```

CI runs these checks on Python 3.11. They cover configuration validation,
rendering, mesh compatibility, failure recovery, and run orchestration. Native
mesh conversion and the coupled solve require the Docker image.

`channel` defines the rectangular flow channel and is required:

- `length`, `height`, and `width` must be positive.
- `wall_thickness` is the positive thickness used for the electrode and
  insulating shell geometry.

`mesh` is required. `target_element_size` is a positive length in metres that
sets the uniform target characteristic length for the Gmsh tetrahedral mesh.
The resulting edge lengths can vary around this target as Gmsh preserves the
geometry and element quality. The same generated mesh is converted for both
OpenFOAM and Elmer.

Internally, the generator installs this value as a constant Gmsh background
size field. Curvature sizing, geometry-point sizing, and propagation of
boundary sizes into the volume are disabled so that they do not compete with
the configured value. This field is a desired local length, not a strict upper
bound on every tetrahedron edge.

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

`openfoam.files` is required in version 2 and maps a case-relative path to a
file description. Paths must name a file directly below `0/`, `constant/`, or
`system/`. Each file accepts an optional `class` (default: `dictionary`) and an
ordered `entries` mapping. Entries render recursively using OpenFOAM syntax:

- mappings become sub-dictionaries;
- arrays become parenthesized lists;
- booleans become `yes` or `no`;
- strings are emitted as OpenFOAM token expressions;
- `null` emits a valueless directive, such as `$p_rgh;` for dictionary reuse.

This generic representation covers transport properties, run control,
discretization schemes, solver settings, and additional OpenFOAM files without
requiring a renderer change for every new keyword. See the complete checked-in
experiment for the baseline `transportProperties`, `controlDict`, `fvSchemes`,
and `fvSolution` settings.

Unknown keys and invalid values fail during validation instead of being
silently ignored.

Duplicate YAML keys and non-finite numeric values are rejected. Arrays can
contain nested arrays and dictionaries; `null` is only valid as a dictionary
entry. Numeric dictionary values retain their full floating-point precision.

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
2. Renders `case.sif`, all generated OpenFOAM initial fields, and configured
   OpenFOAM files.
3. Copies the maintained base `constant/` and `system/` inputs before replacing
   configuration-owned files with their rendered versions.
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

The provided `linear-hall-coarse.yaml` uses a larger target element size and
also requires `--mesh`:

```sh
mhd prepare /experiments/linear-hall-coarse.yaml --mesh
```

Useful options are:

- `--name NAME` chooses the directory name beneath `/runs`.
- `--ranks N` partitions each solver into `N` MPI ranks; the minimum and default
  are `2`. OpenFOAM's coupled parallel execution rejects a single rank.
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

The checked-in examples use `endTime: 1.0e-8` for short pipeline smoke tests.
Choose the run duration and numerical settings needed for your experiment
before drawing physical conclusions from its output.

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
