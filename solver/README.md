# FVC2 Solver

This directory contains the FVC2 flow solver. The solver reads its inputs from
whatever directory it is invoked from (current working directory). The build
produces a named executable: `fvc2_solver`.

## What the solver is

- 2D compressible **Euler** equations (inviscid, no physical viscosity).
- Finite‑volume method on unstructured triangle/quad meshes.
- Explicit time stepping (forward Euler implemented).
- Ideal‑gas thermodynamics.

Notes:
- `diffusion_coeff` is currently read but **not used** (no viscous terms).
- Wall boundaries are **slip walls** (no‑penetration; no no‑slip model).
- `time_discretization = 2` (RK4) is reserved but currently not implemented.
- `subsonic_outlet(...)` boundary condition is reserved but currently not implemented.

## Run directory layout

A typical run directory should contain:

- `case.fvc` (mesh + physical names from `msh2fvc`)
- `case.vars` (flow solver control data)

Outputs are written to `output/` in the run directory:
- Tecplot ASCII: `output/soln_c*.dat`
- VTK legacy: `output/soln_c*.vtk`
- VTK XML (when `animate_output: on`): `output/soln_c*.vtu` + `output/soln.pvd`

The solver creates `output/` if it does not exist.

## Build

```sh
make
```

The executable is `fvc2_solver`.

## case.vars format

The `case.vars` file uses **key:value** pairs. Lines may include comments starting
with `!` or `#` (inline or full-line). Keys are case‑insensitive; values are read
with Fortran list‑directed input. Unknown keys are an error.

### Required keys

- `time_step`: time step size
- `output_interval`: output cadence (steps)
- `max_steps`: number of steps to run
- `gamma`: ratio of specific heats
- `molar_mass`: molar mass (used to set `rspec = R / molar_mass`)
- `diffusion_coeff`: diffusion coefficient (stored as `k`, **not used yet**)
- `time_discretization`: `1` = forward Euler, `2` = RK4 (reserved; currently not implemented)
- `convective_scheme`: `0` = off, `1` = AUSM+, `2` = Roe
- `boundary_conditions`: a block that maps **mesh physical line names** to
  boundary condition types and optional values.

The solver matches these names against the **1D physical names** in the
`[phys_names]` section of `case.fvc` (written by `msh2fvc`). All **1D physical
names** must appear in the block.

Input validation:
- `time_step > 0`
- `output_interval > 0`
- `max_steps >= 0`
- `gamma > 1`
- `molar_mass > 0`
- `diffusion_coeff >= 0`
- `time_discretization` must be `1` or `2`
- `convective_scheme` must be `0`, `1`, or `2`

Boundary condition types:
- `wall()` — slip wall (no‑penetration)
- `inlet(rho p u v)` — prescribed primitive inflow state
- `extrapolated()` — zero‑gradient/outflow
- `farfield(rho p u v)` — inflow/outflow based on normal velocity; inflow uses
  the provided primitive state
- `subsonic_outlet(p)` — reserved keyword (currently not implemented)

Boundary value packing in `bc_values`:
- `inlet(rho p u v)` and `farfield(rho p u v)` map to
  `bc_values(1:4,:) = [rho, p, u, v]`
- `subsonic_outlet(p)` maps to `bc_values(1,:) = p` (reserved; currently not implemented)

Example (all 1D physical names must be listed):

```
boundary_conditions:<
inlet=inlet(1.0 1.0 5.856 0.0)
outlet=extrapolated()
walls=wall()
>
```

### Optional initial-condition keys

- `initial_condition_type`: `uniform`, `fwd_step`, or `shock_tube` (also accepts `1/2/3`)
- `ic_uniform_state`: `rho p u v` (required if `initial_condition_type: uniform`)
- `shock_split_x`: interface x-location for shock tube
- `shock_left_state`:  `rho p u v`
- `shock_right_state`: `rho p u v`

Defaults (if not overridden):
- `initial_condition_type = uniform`
- `ic_uniform_state = 1.4 1.0 3.0 0.0`
- `shock_split_x = 0.5`
- `shock_left_state = 1.0 1.0 0.0 0.0`
- `shock_right_state = 0.125 0.1 0.0 0.0`

### Optional output format

- `output_format_mode`: `tecplot`, `vtk`, or `both` (also accepts `0/1/2`)
- `animate_output`: `on/off`, `true/false`, `yes/no`, `1/0`

If omitted, the defaults are `output_format_mode = both` and `animate_output = off`.
When `animate_output` is on and VTK output is enabled, the solver writes `.vtu`
files and a `soln.pvd` time series. Tecplot output is still written when
`output_format_mode` is `tecplot` or `both`.
