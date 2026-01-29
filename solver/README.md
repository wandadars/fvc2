# FVC2 Solver

This directory contains the FVC2 flow solver. The solver reads its inputs from
whatever directory it is invoked from (current working directory). The build
produces a named executable: `fvc2_solver`.

## Run directory layout

A typical run directory should contain:

- `prj.ini` (run control file)
- `grid/` (mesh files from `sort_grid`: `nodel`, `elnorm`, `elml`, `edgl`)

If `grid/` is not present, the solver falls back to reading the mesh files from the
current working directory.

Outputs are written to `output/` in the run directory:
- Tecplot ASCII: `output/soln_c*.dat`
- VTK legacy: `output/soln_c*.vtk`

The solver creates `output/` if it does not exist.

## Build

```sh
make
```

The executable is `fvc2_solver`.

## prj.ini format

The `prj.ini` file is read in fixed order. Lines may include comments starting
with `!` or `#` (inline or full-line). Values are read line-by-line in this
sequence:

1) `dt`
2) `iostep`
3) `nstepmax`
4) `rgam`
5) `mmass`
6) `k`
7) `time_scheme`
8) `conv_scheme`
9) `bcvals(1..4,1)`  (Dirichlet)
10) `bcvals(1..4,2)`
11) `bcvals(1..4,3)`
12) `bcvals(1..4,4)`
13) `bcvals(1..4,5)`
14) `bcvals(1..2,6)`  (Flux)
15) `bcvals(1..2,7)`
16) `bcvals(1..2,8)`
17) `bcvals(1..2,9)`
18) `bcvals(1..2,10)`

For supersonic inflow boundaries (type `2` from the mesh tags), the solver uses
`bcvals(1..4,2)` as **primitive** inflow state: `rho p u v`. If these are not set
or are non‑physical, it falls back to `ic_uniform`.

### Optional initial-condition block

If present, the following 5 lines must appear **after** the BC blocks:

19) `ic_type`  (1=uniform, 2=fwd_step legacy, 3=shock_tube)
20) `ic_uniform`: `rho p u v`
21) `ic_split_x` (interface x-location for shock tube)
22) `ic_left`:  `rho p u v`
23) `ic_right`: `rho p u v`

If the IC block is omitted entirely, defaults are used (uniform with
`rho=1.4, p=1.0, u=3.0, v=0.0`, and shock-tube defaults for type 3).

### Optional output format

If present, the next line selects output format:

24) `output_mode` (0=Tecplot only, 1=VTK only, 2=both)

If omitted, the default is `2` (both).
