#!/usr/bin/env julia
#
# Quick-look plots of FastIsostasy's basic (non-extended) test output for one
# output version.
#
# Usage:
#   julia --project=. scripts/plot_version.jl <version> [experiment...]
#
# Examples:
#   julia --project=. scripts/plot_version.jl v0.3
#       -> plots every bedtest_<experiment>.nc found under output/v0.3
#   julia --project=. scripts/plot_version.jl v0.3 test1a test4b
#       -> plots only those two experiments
#
# Each figure (heatmaps of H_ice, w_viscous, dwdt, z_bed at the last saved
# time) is written to plot-tests/plots/<version>/<experiment>.png.

include(joinpath(@__DIR__, "..", "src", "PlotTests.jl"))
using .PlotTests

function main(args)
    if isempty(args)
        println(stderr, "Usage: julia --project=. scripts/plot_version.jl <version> [experiment...]")
        println(stderr, "Available versions: ", join(list_versions(), ", "))
        return 1
    end

    version = args[1]
    experiments = length(args) > 1 ? args[2:end] : list_experiments(version)

    println("Plotting $(length(experiments)) experiment(s) from output/$(version)...")
    plot_version(version; experiments)
    println("Done.")
    return 0
end

exit(main(ARGS))
