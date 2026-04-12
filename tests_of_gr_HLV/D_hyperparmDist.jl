# Usage:
#   julia --project=. D_hyperparmDist.jl <config_file.toml> [N] [waveform_family]
#
# Arguments:
#   <config_file.toml> : Required analysis configuration file.
#   [N]                : Optional number of events to draw. Defaults to `n_catalog`
#                        from the config file.
#   [waveform_family]  : Optional waveform family to use. Defaults to the first
#                        waveform family listed in the config file.
#
# The script loads a random subset of events from the Fisher-results file that:
#   1. pass the configured `snr` threshold,
#   2. pass the configured `isnr` threshold,
#   3. have an invertible Fisher matrix,
# for the chosen waveform family and across all configured PN orders.
#
# The random selection is seeded with `config["seed"]`, and the script returns
# the selected `delta_k` and `dphi_k` values for each configured PN order.

using TOML
using HDF5
using Random
using CairoMakie
using Printf

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

const GRID_LIMS = Dict(
    "-1"       => (-0.01,0.01,0.01),
    "0"        => (-0.5,0.5,0.5),
    "0.5"      => (-1.0,1.0,1.0),
    "1"        => (-0.8,0.8,0.8),
    "1.5"      => (-0.5,0.5,0.5),
    "2"        => (-5.0,5.0,5.0),
    "log(2.5)" => (-2.0,2.0,2.0),
    "3"        => (-2.0,2.0,2.0),
    "log(3.)"  => (-15.0,15.0,15.0),
    "3.5"      => (-8.0,8.0,8.0),
)
const CI_LEVELS_2D = (
    1.0 - exp(-0.5),
    1.0 - exp(-2.0),
    1.0 - exp(-4.5),
)
const CONTOUR_FILL_COLORS = [:aliceblue, :lightblue, :cornflowerblue]
const CONTOUR_LINE_COLOR = :royalblue4
const GR_CONTOUR_COLOR = :black
const HYPERPARAM_FIGURE_SIZE = (900, 700)
const HYPERPARAM_PANEL_GAP = 10
const GR_REFERENCE_COLOR = :darkorange

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_hyperparam_plot_output_files(config::Dict, pno::AbstractString)
    output_dir = joinpath(
        @__DIR__,
        config["plot_outdir"],
        "distributions",
        createSED.pnoString(pno),
    )

    return (
        joinpath(output_dir, "hyperparam_dist.png"),
        joinpath(output_dir, "hyperparam_dist.pdf"),
    )
end

"""
Build a `1000 x 1000` `(mu, sigma)` grid from the limits specified in
`GRID_LIMS[pno]`.

The tuple is interpreted as `(mu_min, mu_max, sigma_max)`. The `sigma` grid is
always bounded from below by `0`.
"""
function build_hyperparam_grid(pno::AbstractString)
    haskey(GRID_LIMS, pno) ||
        throw(ArgumentError("Missing grid limits for PN order `$(pno)` in `GRID_LIMS`."))

    mu_min, mu_max, sigma_max = Float64.(GRID_LIMS[pno])
    mu_min < mu_max || throw(ArgumentError("For PN order `$(pno)`, `mu_min` must be smaller than `mu_max`."))
    sigma_max > 0.0 || throw(ArgumentError("For PN order `$(pno)`, `sigma_max` must be positive."))

    mu_spacing = (mu_max - mu_min) / 999.0
    mu_neg = collect(0.0:-mu_spacing:mu_min)
    mu_pos = collect(0.0:mu_spacing:mu_max)
    mu = vcat(reverse(mu_neg[2:end]), mu_pos)
    sigma = collect(range(0.0, sigma_max; length=1000))

    0.0 in sigma || throw(ArgumentError("For PN order `$(pno)`, the generated sigma-grid does not contain 0.0 exactly. Adjust `GRID_LIMS` so that 0 lies on the grid."))

    return mu, sigma
end

"""
Load a seeded random subset of hyperparameter-summary data from the Fisher
results file for one waveform family.

The function first identifies the events that satisfy all configured selection
requirements across every PN order listed in `config["pn_orders"]`:
`snr > config["snr_threshold"]`, `isnr > config["snr_inspiral_threshold"]`,
and an invertible Fisher matrix. It then draws `n_catalog` events uniformly at
random from that eligible set using `config["seed"]`.

For each configured PN order, the function returns the selected `delta_k` and
`dphi_k` values corresponding to the sampled event indices.
"""
function load_selected_hyperparam_data(
    fisher_results_file::AbstractString,
    config::Dict,
    n_catalog::Int,
    waveform_family::AbstractString,
    )

    waveform_family in config["waveform_families"] ||
        throw(ArgumentError("Unknown waveform family `$(waveform_family)`. Expected one of $(config["waveform_families"])."))
    n_catalog > 0 || throw(ArgumentError("`n_catalog` must be positive."))

    rng = MersenneTwister(config["seed"])
    selected_indices = Int[]
    selected_delta_k = Dict{String, Vector{Float64}}()
    selected_dphi_k = Dict{String, Vector{Float64}}()

    h5open(fisher_results_file, "r") do file
        eligible_index = nothing

        for pno in config["pn_orders"]
            pno_group = createSED.pnoString(pno)
            wf_group = file[pno_group][config["network"]][waveform_family]

            snr = read(wf_group, "snr")
            isnr = read(wf_group, "isnr")
            invc = read(wf_group, "invc")

            selection_index = BitVector(
                (snr .> config["snr_threshold"]) .&
                (isnr .> config["snr_inspiral_threshold"]) .&
                invc
            )

            if isnothing(eligible_index)
                eligible_index = selection_index
            else
                eligible_index .&= selection_index
            end
        end

        n_eligible = sum(eligible_index)
        n_catalog <= n_eligible ||
            throw(ArgumentError("Requested N=$(n_catalog) events, but only $(n_eligible) events pass the cuts for waveform `$(waveform_family)`."))

        eligible_event_indices = findall(eligible_index)
        sampled_positions = randperm(rng, n_eligible)[1:n_catalog]
        selected_indices = sort(eligible_event_indices[sampled_positions])

        for pno in config["pn_orders"]
            pno_group = createSED.pnoString(pno)
            wf_group = file[pno_group][config["network"]][waveform_family]

            delta_k = read(wf_group, "delta_k")
            dphi_k = read(wf_group, "dphi_k")

            selected_delta_k[pno] = delta_k[selected_indices]
            selected_dphi_k[pno] = dphi_k[selected_indices]
        end
    end

    return selected_indices, selected_delta_k, selected_dphi_k
end

function save_hyperparam_distribution_plot(
    config::Dict,
    pno::AbstractString,
    mu_grid::Vector{Float64},
    sigma_grid::Vector{Float64},
    p_mu_sigma::Matrix{Float64},
    p_sigma::Vector{Float64},
    p_mu::Vector{Float64},
    contour_info::Dict{Float64, NamedTuple{(:level, :mu_lines, :sigma_lines), Tuple{Float64, Vector{Vector{Float64}}, Vector{Vector{Float64}}}}},
    gr_level::Float64,
    summary_text::AbstractString,
    )

    contour_levels = sort([contour_info[CI].level for CI in CI_LEVELS_2D])
    fill_levels = vcat(contour_levels, [maximum(p_mu_sigma) + eps()])

    fig = Figure(size=HYPERPARAM_FIGURE_SIZE, backgroundcolor=:white)
    top_ax = Axis(
        fig[1, 1],
        backgroundcolor=:white,
        xlabel="",
        ylabel="p(μ|D)",
        title="PN order: $(pno)",
    )
    main_ax = Axis(
        fig[2, 1],
        backgroundcolor=:white,
        xlabel="μ",
        ylabel="σ",
    )
    side_ax = Axis(
        fig[2, 2],
        backgroundcolor=:white,
        xlabel="p(σ|D)",
        ylabel="",
    )
    info_ax = Axis(
        fig[1, 2],
        backgroundcolor=:white,
    )

    contourf!(main_ax, mu_grid, sigma_grid, p_mu_sigma; levels=fill_levels, colormap=CONTOUR_FILL_COLORS)
    contour!(main_ax, mu_grid, sigma_grid, p_mu_sigma; levels=contour_levels, color=CONTOUR_LINE_COLOR, linewidth=1.5)
    contour!(main_ax, mu_grid, sigma_grid, p_mu_sigma; levels=[gr_level], color=GR_CONTOUR_COLOR, linewidth=1.8)

    lines!(top_ax, mu_grid, p_mu; color=CONTOUR_LINE_COLOR, linewidth=2.0)
    lines!(side_ax, p_sigma, sigma_grid; color=CONTOUR_LINE_COLOR, linewidth=2.0)
    vlines!(main_ax, [0.0]; color=GR_REFERENCE_COLOR, linewidth=2.0)
    vlines!(top_ax, [0.0]; color=GR_REFERENCE_COLOR, linewidth=2.0)
    text!(info_ax, 0.03, 0.97; text=summary_text, space=:relative, align=(:left, :top), color=:black)

    xlims!(main_ax, minimum(mu_grid), maximum(mu_grid))
    ylims!(main_ax, minimum(sigma_grid), maximum(sigma_grid))
    xlims!(top_ax, minimum(mu_grid), maximum(mu_grid))
    ylims!(top_ax, 0.0, 1.05 * maximum(p_mu))
    xlims!(side_ax, 0.0, 1.05 * maximum(p_sigma))
    ylims!(side_ax, minimum(sigma_grid), maximum(sigma_grid))
    xlims!(info_ax, 0.0, 1.0)
    ylims!(info_ax, 0.0, 1.0)

    hidexdecorations!(top_ax, grid=false)
    hideydecorations!(side_ax, grid=false)
    hidedecorations!(info_ax)
    hidespines!(info_ax)

    rowsize!(fig.layout, 1, Relative(0.25))
    rowsize!(fig.layout, 2, Relative(0.75))
    colsize!(fig.layout, 1, Relative(0.78))
    colsize!(fig.layout, 2, Relative(0.22))
    rowgap!(fig.layout, HYPERPARAM_PANEL_GAP)
    colgap!(fig.layout, HYPERPARAM_PANEL_GAP)

    png_output_file, pdf_output_file = get_hyperparam_plot_output_files(config, pno)
    mkpath(dirname(png_output_file))

    save(png_output_file, fig)
    save(pdf_output_file, fig)

    return png_output_file, pdf_output_file
end

function run_hyperparam_dist(
    config::Dict;
    n_catalog::Int=config["n_catalog"],
    waveform_family::String=first(config["waveform_families"]),
    )
    fisher_results_file = get_fisher_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    println("Using Fisher results file: $(fisher_results_file)")
    println("Waveform family: $(waveform_family)")
    println("Requested number of events: $(n_catalog)")
    println("Using random seed: $(config["seed"])")

    idx_selected, delta_k_dict, dphi_k_dict = load_selected_hyperparam_data(
        fisher_results_file,
        config,
        n_catalog,
        waveform_family,
    )

    println("Selected $(length(idx_selected)) events after applying the SNR and invertibility cuts.")

    contour_dict = Dict{String, Dict{Float64, NamedTuple{(:level, :mu_lines, :sigma_lines), Tuple{Float64, Vector{Vector{Float64}}, Vector{Vector{Float64}}}}}}()
    for pno in config["pn_orders"] #["0"]#
        hyperparam_dist = HierDist.hyperparamDistTIGER(
            dphi_k_dict[pno],
            delta_k_dict[pno],
        )
        mu_grid, sigma_grid = build_hyperparam_grid(pno)

        p_mu_sigma, p_sigma, p_mu, _, _ = HierDist.getDistributionOnGrid(
            mu_grid,
            sigma_grid,
            hyperparam_dist,
        )

        idx_mu_zero = findfirst(==(0.0), mu_grid)
        idx_sigma_zero = findfirst(==(0.0), sigma_grid)
        isnothing(idx_mu_zero) && throw(ErrorException("The mu grid does not contain 0.0."))
        isnothing(idx_sigma_zero) && throw(ErrorException("The sigma grid does not contain 0.0."))

        gr_level = p_mu_sigma[idx_mu_zero, idx_sigma_zero]
        gr_quantile = HierDist.getEnclosedIsoVolProbOGD(p_mu_sigma, gr_level)

        mu_q05, mu_q50, mu_q95 = HierDist.quantile1dOGD(mu_grid, p_mu, [0.05, 0.5, 0.95])
        sigma_q90 = HierDist.quantile1dOGD(sigma_grid, p_sigma, 0.90)
        summary_text = @sprintf(
            "μ = %.3f +%.3f / -%.3f (90%% CI)\nσ < %.3f (90%%)\nGR_quantile = %.2f%%%%",
            mu_q50,
            mu_q95 - mu_q50,
            mu_q50 - mu_q05,
            sigma_q90,
            100.0 * gr_quantile,
        )

        contour_dict[pno] = Dict{Float64, NamedTuple{(:level, :mu_lines, :sigma_lines), Tuple{Float64, Vector{Vector{Float64}}, Vector{Vector{Float64}}}}}()
        for CI in CI_LEVELS_2D
            contour_level = HierDist.getContourLevelOGD(p_mu_sigma, CI)
            contour_lines_mu, contour_lines_sigma = HierDist.getCredibleContourOGD(
                mu_grid,
                sigma_grid,
                p_mu_sigma,
                CI,
            )
            contour_dict[pno][CI] = (
                level=contour_level,
                mu_lines=contour_lines_mu,
                sigma_lines=contour_lines_sigma,
            )
        end

        png_output_file, pdf_output_file = save_hyperparam_distribution_plot(
            config,
            pno,
            mu_grid,
            sigma_grid,
            p_mu_sigma,
            p_sigma,
            p_mu,
            contour_dict[pno],
            gr_level,
            summary_text,
        )

        println("  PN order $(pno): built hyperparamDistTIGER, evaluated the distribution on a 1000x1000 mu-sigma grid, calculated the 1/2/3 sigma 2D contours, and saved")
        println("    $(png_output_file)")
        println("    $(pdf_output_file)")
    end

    return (
        selected_indices=idx_selected,
        selected_delta_k=delta_k_dict,
        selected_dphi_k=dphi_k_dict,
        contours_2d=contour_dict,
    )
end

#----------------------------------------------------------------------------#
# RUN MAIN
#----------------------------------------------------------------------------#
function main(args=ARGS)
    if !(1 <= length(args) <= 3)
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml> [N] [waveform_family]"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    config = read_config(config_file)

    n_catalog = length(args) >= 2 ? parse(Int, args[2]) : config["n_catalog"]
    waveform_family = length(args) >= 3 ? String(args[3]) : first(config["waveform_families"])

    return run_hyperparam_dist(
        config;
        n_catalog=n_catalog,
        waveform_family=waveform_family,
    )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
