# FVC2

This repo contains:
- `msh2fvc/` mesh converter (Gmsh `.msh` → `case.fvc`)
- `solver/` 2D finite‑volume solver
- `test_cases/` example cases

## Build (solver)

```sh
make -C solver        # optimized build
make -C solver debug  # debug build with runtime checks
```

## Profiling (solver)

Build with gprof instrumentation:

```sh
make -C solver profile
```

Run a case (gprof writes `gmon.out` in the run directory):

```sh
cd test_cases/sod
../../solver/fvc2_solver
```

Generate a report:

```sh
gprof ../../solver/fvc2_solver gmon.out > gprof.txt
```

### Visualizing gprof output

Two common options:

1) **Call graph (SVG) via gprof2dot**

```sh
# install: pip install gprof2dot
gprof2dot -f gprof gprof.txt | dot -Tsvg -o gprof.svg
# then open gprof.svg in a browser
```

2) **KCachegrind / QCachegrind**

```sh
# install: sudo apt install kcachegrind   (or qcachegrind)
# convert:
gprof2calltree gprof.txt > gprof.calltree
# view:
kcachegrind gprof.calltree
```

Notes:
- `gprof2dot` requires Graphviz (`dot`). Install via `sudo apt install graphviz`.
- If you used `GMON_OUT_PREFIX`, pass the specific `gmon` file to `gprof`.

Optional: keep profiles per run using `GMON_OUT_PREFIX`:

```sh
GMON_OUT_PREFIX=output/gmon ../../solver/fvc2_solver
```

If you prefer Linux `perf`:

```sh
perf record ../../solver/fvc2_solver
perf report
```
