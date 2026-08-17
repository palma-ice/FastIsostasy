module PlotTests

using NCDatasets
using CairoMakie

export list_versions, list_experiments, plot_experiment, plot_version

# plot-tests/src -> plot-tests -> FastIsostasy
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const OUTPUT_ROOT = joinpath(REPO_ROOT, "output")
const PLOTS_ROOT = normpath(joinpath(@__DIR__, "..", "plots"))

# Fields shown by default: the primary isostasy outputs in the *basic*
# bedtest_<experiment>.nc files (not the *_extended ones, which carry a much
# larger set of diagnostics and are intentionally skipped here for speed).
const DEFAULT_FIELDS = (:H_ice, :w_viscous, :dwdt, :z_bed)

"""
    list_versions() -> Vector{String}

List output version folders directly under `FastIsostasy/output` (e.g.
"v0.1", "v0.2", "v0.3"). Non-version entries (such as "not-released") are
included as-is; this just reflects whatever directories exist.
"""
list_versions() = sort(filter(isdir ∘ Base.Fix1(joinpath, OUTPUT_ROOT), readdir(OUTPUT_ROOT)))

"""
    list_experiments(version) -> Vector{String}

List experiment codes (e.g. "test1a") for which a *basic* (non-extended)
`bedtest_<experiment>.nc` file exists under `output/<version>`.
"""
function list_experiments(version::AbstractString)
    dir = joinpath(OUTPUT_ROOT, version)
    isdir(dir) || error("No such output directory: $dir")
    files = filter(readdir(dir)) do f
        startswith(f, "bedtest_") && endswith(f, ".nc") && !endswith(f, "_extended.nc")
    end
    isempty(files) && error("No basic bedtest_*.nc files found in $dir")
    experiments = [f[length("bedtest_")+1:end-length(".nc")] for f in files]
    return sort(experiments)
end

"""
    plot_experiment(version, experiment; fields = $(DEFAULT_FIELDS), itime = :last)

Quick-look figure for one experiment's basic output file
(`output/<version>/bedtest_<experiment>.nc`): a heatmap of each of `fields`
at time index `itime` (`:last` for the final saved time, or an Int), each
with its own colorbar. Saved to `plot-tests/plots/<version>/<experiment>.png`
and the path is returned.
"""
function plot_experiment(version::AbstractString, experiment::AbstractString;
        fields = DEFAULT_FIELDS, itime = :last)

    file = joinpath(OUTPUT_ROOT, version, "bedtest_$(experiment).nc")
    isfile(file) || error("No such file: $file")

    # Fixed pixel sizes (rather than relying on GridLayout to shrink columns
    # around a DataAspect()-locked axis, which it won't do on its own) so the
    # layout -- and resize_to_layout! below -- has no slack to fill.
    panel = 260

    # Start from a minimal canvas: by default GridLayout's Auto() columns/rows
    # stretch to fill any leftover figure area, which is exactly the
    # whitespace we're trying to avoid. Starting small leaves nothing to
    # stretch into, so the fixed panel/colorbar sizes above are what determine
    # the final geometry once resize_to_layout! runs.
    fig = Figure(size = (10, 10), figure_padding = 6)

    NCDataset(file) do ds
        xc = ds["xc"][:]
        yc = ds["yc"][:]
        t  = ds["time"][:]
        it = itime === :last ? lastindex(t) : itime

        Label(fig[0, 1:2*length(fields)],
            "$(experiment)  ($(version))   ·   t = $(round(t[it], digits=1)) yr";
            fontsize = 16)

        for (i, var) in enumerate(fields)
            row, col = fldmod1(i, 2)
            data = Array(ds[String(var)][:, :, it])
            ax = Axis(fig[row, 2col-1];
                title = String(var), xlabel = "x (km)", ylabel = "y (km)",
                aspect = DataAspect(), width = panel, height = panel)
            hm = heatmap!(ax, xc, yc, data)
            Colorbar(fig[row, 2col], hm; width = 10, height = panel)
        end
    end

    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 8)
    resize_to_layout!(fig)

    outdir = joinpath(PLOTS_ROOT, version)
    mkpath(outdir)
    outfile = joinpath(outdir, "$(experiment).png")
    save(outfile, fig)
    return outfile
end

"""
    plot_version(version; fields = $(DEFAULT_FIELDS), experiments = list_experiments(version))

Run [`plot_experiment`](@ref) for every basic output file found under
`output/<version>` (or just `experiments`, if given), printing progress and
returning the vector of saved file paths.
"""
function plot_version(version::AbstractString; fields = DEFAULT_FIELDS,
        experiments = list_experiments(version))

    outfiles = String[]
    for experiment in experiments
        outfile = plot_experiment(version, experiment; fields)
        println("  wrote ", outfile)
        push!(outfiles, outfile)
    end
    return outfiles
end

end # module PlotTests
