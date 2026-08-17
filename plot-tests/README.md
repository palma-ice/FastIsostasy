# plot-tests

Quick-look plots of FastIsostasy's test output (`FastIsostasy/output/<version>/bedtest_*.nc`),
using CairoMakie + NCDatasets. Only the basic output files are plotted, not
the `*_extended.nc` ones.

## Usage

```sh
cd FastIsostasy/plot-tests
julia --project=. scripts/plot_version.jl v0.3            # all experiments in output/v0.3
julia --project=. scripts/plot_version.jl v0.3 test1a     # just one
```

Figures are written to `plots/<version>/<experiment>.png`: a 2x2 grid of
`H_ice`, `w_viscous`, `dwdt`, `z_bed` at the last saved time.

Reusable functions (`list_versions`, `list_experiments`, `plot_experiment`,
`plot_version`) live in `src/PlotTests.jl`.

## Cluster note: `LD_LIBRARY_PATH`

This cluster's module system (netcdf, cairo, glib, openssl, ...) sets
`LD_LIBRARY_PATH`, which shadows Julia's own bundled artifact libraries
(`libssl.so`, `libgobject-2.0.so`, ...) with incompatible system versions and
breaks precompilation/loading of NCDatasets and CairoMakie. Run `Pkg` and
`julia` commands here with it cleared:

```sh
env -u LD_LIBRARY_PATH julia --project=. scripts/plot_version.jl v0.3
```
