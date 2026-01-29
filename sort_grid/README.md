# FVC2 Mesh Utility (sort_grid)

This directory contains the mesh pre-processor used by FVC2. It converts a Gmsh `.msh` file into the connectivity and geometry files expected by the solver. Both legacy MSH 2.x and modern MSH 4.x formats are supported. This utility generates four ASCII output files containing mesh information: `nodel`, `elnorm`, `elml`, and `edgl`.

## Workflow overview

1) Create a 2D geometry in Gmsh (z = 0 for points/lines).
2) Define physical entities:
   - The 2D plane surface must be a physical surface (name not important).
   - Boundary curves must be physical lines so they mesh as line elements.
   - Boundary naming uses the name prefix to pick BC type:
     - wall: `wall`, `bump`, `top_`, `symm` -> wall (type 1)
     - inflow: `infl`, `inle` -> supersonic inflow (type 2)
     - outflow: `outf`, `outl` -> supersonic outflow (type 3)
   - If the physical surface is missing, Gmsh will only write line elements and
     `sort_grid` will produce zero cells.
3) Apply geometric constraints (element type, transfinite lines, etc.).
   - Only triangles and quadrilaterals are recognized by FVC2.
4) Generate the mesh and export a `.msh` file.

## Running the utility

1) Set the mesh filename in `prj.ini`. Optional lines allow control of output:
   - `write_vtk_mesh` (0=skip VTK, 1=write VTK, default 1)
   - `verbosity` (0=quiet, 1=summary, 2=per-element, default 1)
2) Build and run:

```sh
make
./sort_grid
```

3) (Optional) copy outputs into the solver directory:

```sh
make move_all
```

Outputs are ASCII files written in this directory: `nodel`, `elnorm`, `elml`, `edgl`.
A debug VTK file `grid_cXXXXX.vtk` is also written.

## Output file formats

For ease of output, each file is written as a single column vector with implied
striding (i.e., blocks of fixed size).

### Node file: `nodel`

Purpose: node numbering and coordinates.

- First number: `nn` = total number of nodes.
- Then `nn` blocks of two numbers: x, y of each node.
- Natural ordering: block 1 is node 1, block 2 is node 2, etc.

Stored in `nodel(i,j)` with `j = 1..nn`, and `i = 1` (x) or `i = 2` (y).

### Normal/geometry file: `elnorm`

Purpose: element centroids, areas, and outward face normals.

- First number: `nel` = total number of 2D elements.
- Then `nel` blocks of eleven numbers:
  - 1-2: centroid x, y
  - 3: element area
  - 4-5: face 1 normal (x, y)
  - 6-7: face 2 normal (x, y)
  - 8-9: face 3 normal (x, y)
  - 10-11: face 4 normal (x, y) (quads only)

Stored in `elnorm(i,j)` with `j = 1..nel`.

### Element file: `elml`

Purpose: element connectivity.

- First number: `nel`.
- Then `nel` blocks of nine numbers:
  - 1-4: node numbers of the element
  - 5-8: edge numbers of the element (indices into `edgl`)
  - 9: number of sides (3 or 4)

Stored in `elml(i,j)` with `j = 1..nel`.

### Edge file: `edgl`

Purpose: edge connectivity and boundary classification.

- First number: `ned` = total number of edges.
- Then `ned` blocks of six numbers:
  - 1-2: node numbers of the edge
  - 3-4: element numbers on either side (for boundary edges, only 3 is used)
  - 5: boundary type (0 none, 1 wall, 2 supersonic inflow, 3 supersonic outflow)
  - 6: not currently used

Stored in `edgl(i,j)` with `j = 1..ned`.
