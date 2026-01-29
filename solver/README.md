# FVC2 Solver

This directory contains the FVC2 flow solver. The solver reads its inputs from
whatever directory it is invoked from (current working directory). The build
produces a named executable: `fvc2_solver`.

## Run directory layout

A typical run directory should contain:

- `case.fvc` (mesh + physical names from `msh2fvc`)
- `case.vars` (flow solver control data)

Outputs are written to `output/` in the run directory:
- Tecplot ASCII: `output/soln_c*.dat`
- VTK legacy: `output/soln_c*.vtk`

The solver creates `output/` if it does not exist.

## Build

```sh
make
```

The executable is `fvc2_solver`.

## case.vars format

The `case.vars` file uses **key:value** pairs. Lines may include comments starting
with `!` or `#` (inline or full-line). Keys are case-insensitive; values are read
with Fortran list-directed input. Unknown keys are an error.

### Required keys

- `time_step`: time step size
- `iostep`: output cadence (steps)
- `max_steps`: number of steps to run
- `gamma`: ratio of specific heats
- `molar_mass`: molar mass (used to set `rspec = R / molar_mass`)
- `diffusion_coeff`: diffusion coefficient (stored as `k`)
- `time_discretization`: `1` = forward Euler, `2` = RK4 (RK4 not yet implemented)
- `convective_scheme`: `0` = off, `1` = AUSM+, `2` = Roe
- `boundary_conditions`: a block that maps **mesh physical line names** to
  boundary condition types and optional values.

For supersonic inflow boundaries, the solver uses the values provided with the
`inflow`/`inlet`/`dirichlet` type as **primitive** inflow state: `rho p u v`. If
these are not set or are non‑physical, it falls back to `ic_uniform`.

The solver matches these names against the **1D physical names** in the
`[phys_names]` section of `case.fvc` (written by `msh2fvc`). All **1D physical
names** must appear in the block.

Example:

```
boundary_conditions:<
inlet=inflow(1.0 1.0 5.856 0.0)
outlet=extrapolated()
walls=slip()
>
```

### Optional initial-condition keys

- `ic_type`: `uniform`, `fwd_step`, or `shock_tube` (also accepts `1/2/3`)
- `ic_uniform`: `rho p u v` (required if `ic_type: uniform`)
- `ic_split_x`: interface x-location for shock tube
- `ic_left`:  `rho p u v`
- `ic_right`: `rho p u v`

If `ic_type` is omitted, defaults are used (uniform with `rho=1.4, p=1.0, u=3.0, v=0.0`,
and shock-tube defaults for type 3).

### Optional output format

- `output_mode`: `tecplot`, `vtk`, or `both` (also accepts `0/1/2`)

If omitted, the default is `both`.
