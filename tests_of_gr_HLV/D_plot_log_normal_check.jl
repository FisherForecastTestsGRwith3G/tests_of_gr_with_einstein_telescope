# This script plots log-normal consistency checks for the bootstrap constraints
# stored in a population-results HDF5 file.
#
# Usage:
#   julia --project=. D_plot_log_normal_check.jl <config_file.toml>
#
# Expected config entries:
#   - `bootstrap_outdir`: directory containing `population_results_<catalog_tag>.h5`
#   - `plot_outdir`: base directory where the output PDFs will be written
#   - `catalog_tag`, `plot_tag`, `network`, `waveform_families`, `pn_orders`
#
# Input:
#   Reads bootstrap samples from
#   `tests_of_gr_HLV/<bootstrap_outdir>/population_results_<catalog_tag>.h5`.
#
# Output:
#   Writes one PDF per waveform family to
#   `tests_of_gr_HLV/<plot_outdir>/checks/bootstrap/`.
#
# Typical workflow:
#   1. Run `B_population_analysis.jl` to generate the population-results file.
#   2. Run this script with the same TOML configuration file.

using TOML
using HDF5
using CairoMakie
using Statistics
using Distributions
using LaTeXStrings

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const LOG_NORMAL_CHECK_SIZE_PER_PANEL = (320, 420)
const LOG_NORMAL_CHECK_BINS = 30
const LOG_NORMAL_CHECK_HIST_COLOR = :darkorange
const LOG_NORMAL_CHECK_CURVE_COLOR = :black

function get_population_results_file(config::Dict)
    return joinpath(@__DIR__, config["bootstrap_outdir"], "population_results_$(config["bootstrap_tag"]).h5")
end

function get_check_plot_output_file(config::Dict, waveform_family::AbstractString)
    return joinpath(
        @__DIR__,
        config["plot_outdir"],
        "checks",
        "bootstrap",
        "log_normal_check_$(waveform_family)_$(config["plot_tag"]).pdf",
    )
end

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#
"""
    read_bootstrap_constraints(population_results_file, config)

Read bootstrap constraint samples from `population_results_file` for the
configured network, waveform families, and PN orders, returning them as a
nested dictionary indexed by waveform family and PN order.
"""
function read_bootstrap_constraints(population_results_file::AbstractString, config::Dict)
    bootstrap_constraints = Dict{String, Dict{String, Vector{Float64}}}()

    h5open(population_results_file, "r") do file
        network_group = file[config["network"]]
        for wf_fam in config["waveform_families"]
            bootstrap_constraints[wf_fam] = Dict{String, Vector{Float64}}()
            waveform_group = network_group[wf_fam]

            for pno in config["pn_orders"]
                pno_group = waveform_group[createSED.pnoString(pno)]
                bootstrap_constraints[wf_fam][pno] = read(pno_group, "bootstrap_constraints")
            end
        end
    end

    return bootstrap_constraints
end

"""
    plot_bootstrap_checks_for_waveform(config, waveform_family, bootstrap_constraints)

Create and save per-PN-order histogram checks for `waveform_family` by plotting
the `log10` of positive bootstrap constraint samples together with a fitted
normal density.
"""
function plot_bootstrap_checks_for_waveform(config::Dict, waveform_family::String, bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}})
    n_pno = length(config["pn_orders"])
    fig = Figure(
        size=(LOG_NORMAL_CHECK_SIZE_PER_PANEL[1] * n_pno, LOG_NORMAL_CHECK_SIZE_PER_PANEL[2]),
        backgroundcolor=:white,
    )

    for (idx, pno) in enumerate(config["pn_orders"])
        values = bootstrap_constraints[waveform_family][pno]

        positive_values = values[values .> 0.0]
        isempty(positive_values) && throw(ArgumentError("No positive bootstrap constraints found for waveform `$(waveform_family)` and PN order `$(pno)`."))
        log_values = log10.(positive_values)

        # find gaussian approximation
        mu_log = mean(log_values)
        sigma_log = max(std(log_values), sqrt(eps(Float64)))
        x_grid = range(minimum(log_values), maximum(log_values), length=400)
        bin_edges = collect(range(first(x_grid), last(x_grid), length=LOG_NORMAL_CHECK_BINS + 1))

        # calculate JS divergence between distribution and gaussian approximation
        empirical_counts = zeros(Float64, LOG_NORMAL_CHECK_BINS)
        @inbounds for value in log_values
            bin_idx = searchsortedlast(bin_edges, value)
            bin_idx = clamp(bin_idx, 1, LOG_NORMAL_CHECK_BINS)
            if value == bin_edges[end]
                bin_idx = LOG_NORMAL_CHECK_BINS
            end
            empirical_counts[bin_idx] += 1.0
        end
        empirical_probs = empirical_counts ./ length(log_values)

        normal_dist = Normal(mu_log, sigma_log)
        gaussian_probs = zeros(Float64, LOG_NORMAL_CHECK_BINS)
        @inbounds for i in 1:LOG_NORMAL_CHECK_BINS
            gaussian_probs[i] = cdf(normal_dist, bin_edges[i + 1]) - cdf(normal_dist, bin_edges[i])
        end
        gaussian_probs ./= sum(gaussian_probs)

        mixture_probs = 0.5 .* (empirical_probs .+ gaussian_probs)
        js_divergence = 0.5 * (
            sum(ifelse(empirical_probs[i] > 0.0, empirical_probs[i] * log(empirical_probs[i] / mixture_probs[i]), 0.0) for i in eachindex(empirical_probs)) +
            sum(ifelse(gaussian_probs[i] > 0.0, gaussian_probs[i] * log(gaussian_probs[i] / mixture_probs[i]), 0.0) for i in eachindex(gaussian_probs))
        )

        ax = Axis(
            fig[1, idx],
            backgroundcolor=:white,
            xlabel=L"\log_{10}(%$(createSED.pnoLatex(pno)))",
            ylabel=idx == 1 ? "density" : "",
            title="log10(x)",
        )

        hist!(
            ax,
            log_values;
            bins=bin_edges,
            normalization=:pdf,
            color=(LOG_NORMAL_CHECK_HIST_COLOR, 0.75),
            strokecolor=:black,
            strokewidth=1.0,
        )

        lines!(
            ax,
            x_grid,
            pdf.(normal_dist, x_grid);
            color=LOG_NORMAL_CHECK_CURVE_COLOR,
            linewidth=2.5,
        )

        text!(
            ax,
            minimum(log_values) + 0.05 * (maximum(log_values) - minimum(log_values)),
            maximum(pdf.(normal_dist, x_grid)) * 0.92,
            text="JS = $(round(js_divergence; digits=4))",
            align=(:left, :top),
            color=:black,
            fontsize=10,
        )

        if idx != 1
            hideydecorations!(ax; grid=false)
        end
    end

    colgap!(fig.layout, 12)

    output_file = get_check_plot_output_file(config, waveform_family)
    mkpath(dirname(output_file))
    save(output_file, fig)
    println("Saved bootstrap check plot for $(waveform_family) to $(output_file)")
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
function run_plot_check(config::Dict)
    population_results_file = get_population_results_file(config)
    isfile(population_results_file) || throw(ArgumentError("Population results file does not exist: $(population_results_file)"))

    bootstrap_constraints = read_bootstrap_constraints(population_results_file, config)
    for waveform_family in config["waveform_families"]
        plot_bootstrap_checks_for_waveform(config, waveform_family, bootstrap_constraints)
    end

    return bootstrap_constraints
end

function main(args=ARGS)
    if length(args) != 1
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml>"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    config = read_config(config_file)

    return run_plot_check(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    CairoMakie.activate!()
    main()
end
