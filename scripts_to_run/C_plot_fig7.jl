using TOML
using HDF5
using CairoMakie
using LaTeXStrings

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

const FIG7_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
const FIG7_SIZE = (1800, 760)
const GUIDE_FONT_SIZE = 32
const TICK_FONT_SIZE = 32
const LEGEND_FONT_SIZE = 28
const FIG7_TOP_ROW_FRACTION = 0.11
const FIG7_COL_GAP = 20
const FIG7_LABEL_ROW_GAP = 10
const VIOLIN_ALPHA = 0.35
const VIOLIN_EDGE_ALPHA = 0.85

function get_population_results_file(config::Dict)
    return joinpath(
        @__DIR__,
        config["bootstrap_outdir"],
        "population_results_$(config["bootstrap_tag"]).h5"
    )
end

function get_fig7_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"], "fig_7")
    stem = "fig7_$(config["plot_tag"])"
    return (
        joinpath(output_dir, stem * ".png"),
        joinpath(output_dir, stem * ".pdf"),
    )
end

function read_fig7_config(config_file::AbstractString)
    config = read_config(config_file)
    raw_config = TOML.parsefile(abspath(config_file))
    plot_config = get(raw_config, "plots", Dict{String, Any}())
    fig7_config = get(plot_config, "fig7", Dict{String, Any}())

    config["fig7_grid_points"] = Int(get(fig7_config, "grid_points", 700))
    config["fig7_plot_conditioned_distribution"] = Bool(get(fig7_config, "plot_conditioned_distribution", true))
    config["fig7_offset_x_axis"] = Float64(get(fig7_config, "offset_x_axis", 0.16))
    config["fig7_violin_width"] = Float64(get(fig7_config, "violin_width", 0.34))
    config["fig7_grouping"] = haskey(fig7_config, "subplots_pn_order_grouping") ?
        [Int.(group) for group in fig7_config["subplots_pn_order_grouping"]] :
        default_fig7_grouping(length(config["pn_orders"]))
    config["fig7_y_limits"] = haskey(fig7_config, "y_axis_limits") ?
        [Tuple(Float64.(limits)) for limits in fig7_config["y_axis_limits"]] :
        nothing

    return config
end

function default_fig7_grouping(n_pn_orders::Integer)
    n_pn_orders >= 10 && return [[1], collect(2:5), collect(6:n_pn_orders)]
    return [collect(1:n_pn_orders)]
end

function relative_x_offset(index::Integer, total::Integer)
    total == 1 && return 0.0
    return ((index - 1) / (total - 1)) * 2.0 - 1.0
end

function waveform_label(waveform_family::AbstractString)
    waveform_family == "PhenomD" && return L"\text{IMRPhenomD}"
    waveform_family == "PhenomHM" && return L"\text{IMRPhenomHM}"
    return latexstring("\\text{$(waveform_family)}")
end

function top_pno_label(pno::AbstractString)
    if startswith(pno, "log(") && endswith(pno, ")")
        order = chop(chop(pno; head=4); tail=1)
        endswith(order, ".") && (order = chop(order; tail=1))
        return latexstring(order, raw"\,\mathrm{PN}^{(\ell)}")
    end

    return latexstring(pno, raw"\,\mathrm{PN}")
end

function clean_likelihood_data(dphi_k::Vector{Float64}, delta_k::Vector{Float64})
    valid = isfinite.(dphi_k) .& isfinite.(delta_k) .& (delta_k .> 0.0)
    return dphi_k[valid], delta_k[valid]
end

function read_selected_likelihood_data(population_results_file::AbstractString, config::Dict)
    selected_dphi_k = Dict{String, Dict{String, Vector{Float64}}}()
    selected_delta_k = Dict{String, Dict{String, Vector{Float64}}}()

    h5open(population_results_file, "r") do file
        network_group = file[config["network"]]
        for wf_fam in config["waveform_families"]
            selected_dphi_k[wf_fam] = Dict{String, Vector{Float64}}()
            selected_delta_k[wf_fam] = Dict{String, Vector{Float64}}()
            waveform_group = network_group[wf_fam]

            for pno in config["pn_orders"]
                pno_group = waveform_group[createSED.pnoString(pno)]
                dphi, delta = clean_likelihood_data(
                    Float64.(read(pno_group, "selected_dphi_k")),
                    Float64.(read(pno_group, "selected_delta_k")),
                )
                isempty(dphi) &&
                    throw(ArgumentError("No finite selected likelihood data for waveform $(wf_fam), PN order $(pno)."))
                selected_dphi_k[wf_fam][pno] = dphi
                selected_delta_k[wf_fam][pno] = delta
            end
        end
    end

    return selected_dphi_k, selected_delta_k
end

function normalized_density(values::Vector{Float64})
    max_value = maximum(values)
    max_value > 0.0 || throw(ArgumentError("Cannot normalize a zero density profile."))
    return values ./ max_value
end

function posterior_profile(dphi_k::Vector{Float64}, delta_k::Vector{Float64}, grid_points::Integer)
    hyperparam_dist = HierDist.hyperparamDistTIGER(dphi_k, delta_k)
    mu_min, mu_max, sigma_min, sigma_max = HierDist.findOptimalGrid(hyperparam_dist; verbose=false)
    mu_grid = HierDist.buildUniform1dGridWithCenter(mu_min, mu_max, grid_points)
    sigma_grid = HierDist.buildUniform1dGridWithCenter(sigma_min, sigma_max, grid_points)
    _, _, p_mu, _, _ = HierDist.getDistributionOnGrid(mu_grid, sigma_grid, hyperparam_dist)
    conditioned_p_mu = HierDist.getNaiveDistributionOnGrid(mu_grid, hyperparam_dist)

    return (
        mu_grid=mu_grid,
        p_mu=normalized_density(p_mu),
        conditioned_p_mu=normalized_density(conditioned_p_mu),
        q05=HierDist.quantile1dOGD(mu_grid, p_mu, 0.05),
        q95=HierDist.quantile1dOGD(mu_grid, p_mu, 0.95),
    )
end

function build_posterior_profiles(selected_dphi_k, selected_delta_k, config::Dict)
    profiles = Dict{String, Dict{String, NamedTuple}}()
    for wf_fam in config["waveform_families"]
        profiles[wf_fam] = Dict{String, NamedTuple}()
        for pno in config["pn_orders"]
            println("Building posterior profile for $(wf_fam), PN order $(pno)")
            profiles[wf_fam][pno] = posterior_profile(
                selected_dphi_k[wf_fam][pno],
                selected_delta_k[wf_fam][pno],
                config["fig7_grid_points"],
            )
        end
    end
    return profiles
end

function panel_limits(pn_orders::AbstractVector{<:AbstractString}, profiles, waveform_families)
    y_min = Inf
    y_max = -Inf
    for wf_fam in waveform_families
        for pno in pn_orders
            profile = profiles[wf_fam][pno]
            y_min = min(y_min, profile.q05)
            y_max = max(y_max, profile.q95)
        end
    end

    width = y_max - y_min
    width > 0.0 || (width = max(abs(y_max), eps()))
    return (y_min - 0.45 * width, y_max + 0.45 * width)
end

function build_top_label_row!(grid::GridLayout, pn_orders::AbstractVector{<:AbstractString})
    for (idx, pno) in enumerate(pn_orders)
        Label(grid[1, idx], top_pno_label(pno);
            fontsize=TICK_FONT_SIZE,
            tellwidth=false,
            tellheight=false,
            halign=:center,
            valign=:bottom)
        colsize!(grid, idx, Relative(1.0 / length(pn_orders)))
    end
    rowsize!(grid, 1, Relative(1.0))
    return grid
end

function configure_panel_axis!(ax::Axis, pn_orders::AbstractVector{<:AbstractString}, y_limits::Tuple{<:Real, <:Real};
    ylabel="", show_ylabel::Bool=true)
    ax.xlabel = "PN order"
    ax.ylabel = show_ylabel ? ylabel : ""
    ax.xlabelsize = GUIDE_FONT_SIZE
    ax.ylabelsize = GUIDE_FONT_SIZE
    ax.xticks = (collect(1:length(pn_orders)), createSED.pnoLatex.(pn_orders))
    ax.xticklabelsize = TICK_FONT_SIZE
    ax.yticklabelsize = TICK_FONT_SIZE
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = false
    ax.yminorgridvisible = true
    ax.yminorgridcolor = (:gray70, 0.35)
    ax.yminorticks = IntervalsBetween(5)
    ax.xtickalign = 1
    ax.ytickalign = 1
    xlims!(ax, 0.5, length(pn_orders) + 0.5)
    ylims!(ax, y_limits...)
    hlines!(ax, [0.0]; color=(:gray40, 0.8), linestyle=:dash, linewidth=2)
    return ax
end

function add_violin_profile!(ax::Axis, x0::Real, y_grid::Vector{Float64}, density::Vector{Float64}, color;
    max_width::Float64, fill_alpha::Float64=VIOLIN_ALPHA, edge_alpha::Float64=VIOLIN_EDGE_ALPHA)
    widths = max_width .* density
    left_points = [Point2f(x0 - widths[idx], y_grid[idx]) for idx in eachindex(y_grid)]
    right_points = [Point2f(x0 + widths[idx], y_grid[idx]) for idx in reverse(eachindex(y_grid))]
    points = vcat(left_points, right_points)
    poly!(ax, points; color=(color, fill_alpha), strokecolor=(color, edge_alpha), strokewidth=1.5)
    return ax
end

function add_conditioned_outline!(ax::Axis, x0::Real, y_grid::Vector{Float64}, density::Vector{Float64}, color;
    max_width::Float64)
    widths = max_width .* density
    lines!(ax, x0 .- widths, y_grid; color=color, linewidth=2.5)
    lines!(ax, x0 .+ widths, y_grid; color=color, linewidth=2.5)
    return ax
end

function plot_panel!(ax::Axis, pn_orders::AbstractVector{<:AbstractString}, profiles, config::Dict)
    n_waveforms = length(config["waveform_families"])
    for (wf_idx, wf_fam) in enumerate(config["waveform_families"])
        color = FIG7_COLORS[mod1(wf_idx, length(FIG7_COLORS))]
        offset = config["fig7_offset_x_axis"] * relative_x_offset(wf_idx, n_waveforms)

        for (pno_idx, pno) in enumerate(pn_orders)
            profile = profiles[wf_fam][pno]
            x0 = pno_idx + offset
            add_violin_profile!(
                ax,
                x0,
                profile.mu_grid,
                profile.p_mu,
                color;
                max_width=config["fig7_violin_width"],
            )
            if config["fig7_plot_conditioned_distribution"]
                add_conditioned_outline!(
                    ax,
                    x0,
                    profile.mu_grid,
                    profile.conditioned_p_mu,
                    color;
                    max_width=config["fig7_violin_width"],
                )
            end
        end
    end
    return ax
end

function add_fig7_legend!(fig::Figure, target_slot, waveform_families)
    elements = [
        PolyElement(
            color=(FIG7_COLORS[mod1(idx, length(FIG7_COLORS))], VIOLIN_ALPHA),
            strokecolor=(FIG7_COLORS[mod1(idx, length(FIG7_COLORS))], VIOLIN_EDGE_ALPHA),
            strokewidth=1.5,
        )
        for idx in eachindex(waveform_families)
    ]
    labels = waveform_label.(waveform_families)

    Legend(target_slot, elements, labels;
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:bottom,
        margin=(10, 10, 10, 10),
        framevisible=true,
        backgroundcolor=(:white, 0.9),
        labelsize=LEGEND_FONT_SIZE)
    return fig
end

function validate_fig7_grouping(config::Dict)
    n_pn_orders = length(config["pn_orders"])
    grouped_indices = collect(Iterators.flatten(config["fig7_grouping"]))
    all((1 .<= grouped_indices) .& (grouped_indices .<= n_pn_orders)) ||
        throw(ArgumentError("Invalid PN-order index in fig7 subplot grouping."))
    length(unique(grouped_indices)) == length(grouped_indices) ||
        throw(ArgumentError("Duplicate PN-order index in fig7 subplot grouping."))
    return config
end

function build_figure(config::Dict, profiles)
    validate_fig7_grouping(config)
    CairoMakie.activate!()

    grouped_orders = [config["pn_orders"][indices] for indices in config["fig7_grouping"]]
    n_panels = length(grouped_orders)
    fig = Figure(size=FIG7_SIZE, backgroundcolor=:white)

    total_width = sum(length.(grouped_orders))
    for panel_idx in 1:n_panels
        top_grid = GridLayout(fig[1, panel_idx])
        panel_orders = grouped_orders[panel_idx]
        ax = Axis(fig[2, panel_idx], backgroundcolor=:white)
        y_limits = isnothing(config["fig7_y_limits"]) ?
            panel_limits(panel_orders, profiles, config["waveform_families"]) :
            config["fig7_y_limits"][panel_idx]

        colsize!(fig.layout, panel_idx, Relative(length(panel_orders) / total_width))
        build_top_label_row!(top_grid, panel_orders)
        configure_panel_axis!(ax, panel_orders, y_limits;
            ylabel=L"\delta\varphi_p",
            show_ylabel=panel_idx == 1)
        plot_panel!(ax, panel_orders, profiles, config)
    end

    rowsize!(fig.layout, 1, Relative(FIG7_TOP_ROW_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG7_TOP_ROW_FRACTION))
    colgap!(fig.layout, FIG7_COL_GAP)
    rowgap!(fig.layout, FIG7_LABEL_ROW_GAP)
    add_fig7_legend!(fig, fig[2, n_panels], config["waveform_families"])

    return fig
end

function run_plot_fig7(config::Dict)
    population_results_file = get_population_results_file(config)
    isfile(population_results_file) ||
        throw(ArgumentError("Population results file does not exist: $(population_results_file)"))

    selected_dphi_k, selected_delta_k = read_selected_likelihood_data(population_results_file, config)
    profiles = build_posterior_profiles(selected_dphi_k, selected_delta_k, config)
    fig = build_figure(config, profiles)

    png_output_file, pdf_output_file = get_fig7_output_files(config)
    mkpath(dirname(png_output_file))
    save(png_output_file, fig)
    save(pdf_output_file, fig)
    println("Saved figure to $(png_output_file)")
    println("Saved figure to $(pdf_output_file)")

    return png_output_file
end

function main(args=ARGS)
    if length(args) != 1
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml>"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    config = read_fig7_config(config_file)

    return run_plot_fig7(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
