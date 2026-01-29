# FVC2 Mesh Utility (msh2fvc)

This directory contains the mesh pre-processor used by FVC2. It converts a Gmsh
`.msh` file into a `case.fvc` file consumed by the solver. Both legacy MSH 2.x
and modern MSH 4.x formats are supported.

## Workflow overview

1) Create a 2D geometry in Gmsh (z = 0 for points/lines).
2) Define physical entities:
   - The 2D plane surface must be a physical surface (name not important).
   - Boundary curves must be physical lines so they mesh as line elements.
   - Boundary naming must match the names used in `case.vars` within
     `boundary_conditions:`.
   - If the physical surface is missing, Gmsh will only write line elements and
     `msh2fvc` will produce zero cells.
3) Apply geometric constraints (element type, transfinite lines, etc.).
   - Only triangles and quadrilaterals are recognized by FVC2.
4) Generate the mesh and export a `.msh` file.

## Running the utility

1) Build and run from your case directory:

```sh
make
./msh2fvc path/to/case.msh
```
`msh2fvc` only takes the mesh filename. It does **not** read `case.vars`.
Outputs:
- `case.fvc` (mesh + physical names for the solver)
- optional `grid_cXXXXX.vtk` if VTK output is enabled in the code

## case.fvc layout

The file is ASCII with simple section markers:

```
FVC2_CASE 1
[phys_names]
<count>
<dim> <tag> <name>
[end_phys_names]
[mesh]
<nn> <nel> <ned>
nodel
<x y> ... per node
elnorm
<11 values> ... per element
elml
<9 values> ... per element
edgl
<6 values> ... per edge
[end_mesh]
```
