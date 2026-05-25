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

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
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
    config_file = abspath(config_file)
    config = read_config(config_file)
    raw_config = TOML.parsefile(config_file)
    plot_config = get(raw_config, "plots", Dict{String, Any}())
    fig7_config = get(plot_config, "fig7", Dict{String, Any}())

    config["fig7_config_file"] = config_file
    config["fig7_grid_points"] = Int(get(fig7_config, "grid_points", 700))
    config["fig7_plot_conditioned_distribution"] = Bool(get(fig7_config, "plot_conditioned_distribution", true))
    config["fig7_offset_x_axis"] = Float64(get(fig7_config, "offset_x_axis", 0.16))
    config["fig7_violin_width"] = Float64(get(fig7_config, "violin_width", 0.34))
    config["fig7_waveform_family"] = String(get(fig7_config, "waveform_family", first(config["waveform_families"])))
    config["fig7_n_events"] = haskey(fig7_config, "n_events") ? Int(fig7_config["n_events"]) : config["n_events"]
    config["fig7_number_of_events_single_realization"] = haskey(fig7_config, "number_of_events_single_realization") ?
        Int(fig7_config["number_of_events_single_realization"]) :
        nothing
    config["fig7_realization_index"] = Int(get(fig7_config, "realization_index", 1))
    config["fig7_observed_events_direct"] = Bool(get(fig7_config, "n_events_and_number_of_events_single_realization_refer_directly_to_observed_events", false))
    config["fig7_grouping"] = haskey(fig7_config, "subplots_pn_order_grouping") ?
        [Int.(group) for group in fig7_config["subplots_pn_order_grouping"]] :
        default_fig7_grouping(length(config["pn_orders"]))
    config["fig7_y_limits"] = haskey(fig7_config, "y_axis_limits") ?
        [Tuple(Float64.(limits)) for limits in fig7_config["y_axis_limits"]] :
        nothing
    config["fig7_comparison_config_files"] = haskey(fig7_config, "comparison_config_files") ?
        String.(fig7_config["comparison_config_files"]) :
        String[]
    config["fig7_series_labels"] = haskey(fig7_config, "series_labels") ?
        String.(fig7_config["series_labels"]) :
        String[]

    return config
end

function resolve_config_path(path::AbstractString, base_config::Dict)
    isabspath(path) && return path
    return abspath(joinpath(dirname(base_config["fig7_config_file"]), path))
end

function network_series_label(config::Dict)
    network_labels = Dict(
        "ETS" => "ET 10 km",
        "ET_45_15km" => "ET 15 km, 45 deg",
        "ET_0_15km" => "ET 15 km, 0 deg",
        "HLV_O3b" => "HLV O3b",
    )
    return get(network_labels, config["network"], config["network"])
end

function fig7_plot_sources(config::Dict)
    comparison_files = config["fig7_comparison_config_files"]
    labels = config["fig7_series_labels"]
    isempty(comparison_files) && return [
        (
            config=config,
            label=waveform_label(wf_fam),
            waveform_family=wf_fam,
        )
        for wf_fam in config["waveform_families"]
    ]

    isempty(labels) || length(labels) == length(comparison_files) ||
        throw(ArgumentError("`plots.fig7.series_labels` must have the same length as `comparison_config_files`."))

    sources = Vector{NamedTuple}(undef, length(comparison_files))
    for (idx, config_file) in enumerate(comparison_files)
        source_config = read_config(resolve_config_path(config_file, config))
        source_config["fig7_grid_points"] = config["fig7_grid_points"]
        source_config["fig7_n_events"] = config["fig7_n_events"]
        source_config["fig7_number_of_events_single_realization"] = config["fig7_number_of_events_single_realization"]
        source_config["fig7_realization_index"] = config["fig7_realization_index"]
        source_config["fig7_observed_events_direct"] = config["fig7_observed_events_direct"]
        waveform_family = config["fig7_waveform_family"]
        waveform_family in source_config["waveform_families"] ||
            throw(ArgumentError("Waveform family `$(waveform_family)` is not listed in $(config_file)."))
        source_label = isempty(labels) ? network_series_label(source_config) : labels[idx]
        sources[idx] = (
            config=source_config,
            label=latexstring("\\text{$(source_label)}"),
            waveform_family=waveform_family,
        )
    end
    return sources
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

function waveform_summary_selection_from_fisher(file, config::Dict, waveform_family::AbstractString)
    summary_selection = nothing
    for pno in config["pn_orders"]
        wf_group = file[createSED.pnoString(pno)][config["network"]][waveform_family]
        snr = Float64.(read(wf_group, "snr"))
        isnr = Float64.(read(wf_group, "isnr"))
        invc = Bool.(read(wf_group, "invc"))

        pno_selection = BitVector(
            (snr .> config["snr_threshold"]) .&
            (isnr .> config["snr_inspiral_threshold"]) .&
            invc
        )
        summary_selection = isnothing(summary_selection) ? pno_selection : BitVector(summary_selection .& pno_selection)
    end
    return summary_selection
end

function realization_slice(
    full_dphi_k::Vector{Float64},
    full_delta_k::Vector{Float64},
    summary_selection::BitVector,
    config::Dict,
)
    events_per_realization = config["fig7_number_of_events_single_realization"]
    events_per_realization === nothing && return clean_likelihood_data(
        full_dphi_k[summary_selection],
        full_delta_k[summary_selection],
    )

    n_events_requested = min(config["fig7_n_events"], length(full_dphi_k))
    observed_direct = config["fig7_observed_events_direct"]

    if observed_direct
        observed_indices = findall(summary_selection)
        n_events_used = min(n_events_requested, length(observed_indices))
        n_realizations = fld(n_events_used, events_per_realization)
        n_realizations > 0 ||
            throw(ArgumentError("Not enough observed events to build one Figure 7 realization."))
        realization_index = config["fig7_realization_index"]
        1 <= realization_index <= n_realizations ||
            throw(ArgumentError("Requested realization $(realization_index), but only $(n_realizations) are available."))
        first_idx = (realization_index - 1) * events_per_realization + 1
        last_idx = realization_index * events_per_realization
        chosen = observed_indices[first_idx:last_idx]
        println("Using observed-events realization $(realization_index)/$(n_realizations) with $(length(chosen)) observed events.")
        return clean_likelihood_data(full_dphi_k[chosen], full_delta_k[chosen])
    end

    n_realizations = fld(n_events_requested, events_per_realization)
    n_realizations > 0 ||
        throw(ArgumentError("Not enough catalog events to build one Figure 7 realization."))
    realization_index = config["fig7_realization_index"]
    1 <= realization_index <= n_realizations ||
        throw(ArgumentError("Requested realization $(realization_index), but only $(n_realizations) are available."))

    first_idx = (realization_index - 1) * events_per_realization + 1
    last_idx = realization_index * events_per_realization
    chunk_indices = collect(first_idx:last_idx)
    observed_chunk_indices = chunk_indices[summary_selection[chunk_indices]]
    println("Using catalog realization $(realization_index)/$(n_realizations): $(length(observed_chunk_indices)) observed events out of $(events_per_realization) catalog events.")
    return clean_likelihood_data(full_dphi_k[observed_chunk_indices], full_delta_k[observed_chunk_indices])
end

function read_realization_likelihood_data_from_fisher(fisher_results_file::AbstractString, config::Dict, waveform_family::AbstractString)
    selected_dphi_k = Dict{String, Vector{Float64}}()
    selected_delta_k = Dict{String, Vector{Float64}}()

    h5open(fisher_results_file, "r") do file
        summary_selection = waveform_summary_selection_from_fisher(file, config, waveform_family)
        for pno in config["pn_orders"]
            wf_group = file[createSED.pnoString(pno)][config["network"]][waveform_family]
            dphi, delta = realization_slice(
                Float64.(read(wf_group, "dphi_k")),
                Float64.(read(wf_group, "delta_k")),
                summary_selection,
                config,
            )
            isempty(dphi) &&
                throw(ArgumentError("Realization has no finite observed events for waveform $(waveform_family), PN order $(pno)."))
            selected_dphi_k[pno] = dphi
            selected_delta_k[pno] = delta
        end
    end

    return selected_dphi_k, selected_delta_k
end

function clean_likelihood_data(dphi_k::Vector{Float64}, delta_k::Vector{Float64})
    valid = isfinite.(dphi_k) .& isfinite.(delta_k) .& (delta_k .> 0.0)
    return dphi_k[valid], delta_k[valid]
end

function read_selected_likelihood_data(population_results_file::AbstractString, config::Dict, waveform_family::AbstractString)
    selected_dphi_k = Dict{String, Vector{Float64}}()
    selected_delta_k = Dict{String, Vector{Float64}}()

    h5open(population_results_file, "r") do file
        network_group = file[config["network"]]
        waveform_group = network_group[waveform_family]

        for pno in config["pn_orders"]
            pno_group = waveform_group[createSED.pnoString(pno)]
            dphi, delta = clean_likelihood_data(
                Float64.(read(pno_group, "selected_dphi_k")),
                Float64.(read(pno_group, "selected_delta_k")),
            )
            isempty(dphi) &&
                throw(ArgumentError("No finite selected likelihood data for waveform $(waveform_family), PN order $(pno)."))
            selected_dphi_k[pno] = dphi
            selected_delta_k[pno] = delta
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

function build_posterior_profiles(selected_dphi_k, selected_delta_k, source_config::Dict, source_label)
    profiles = Dict{String, NamedTuple}()
    for pno in source_config["pn_orders"]
        println("Building posterior profile for $(source_label), PN order $(pno)")
        profiles[pno] = posterior_profile(
            selected_dphi_k[pno],
            selected_delta_k[pno],
            source_config["fig7_grid_points"],
        )
    end
    return profiles
end

function panel_limits(pn_orders::AbstractVector{<:AbstractString}, profiles, series_ids)
    y_min = Inf
    y_max = -Inf
    for series_id in series_ids
        for pno in pn_orders
            profile = profiles[series_id][pno]
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

function plot_panel!(ax::Axis, pn_orders::AbstractVector{<:AbstractString}, profiles, config::Dict, series_ids)
    n_series = length(series_ids)
    for (series_idx, series_id) in enumerate(series_ids)
        color = FIG7_COLORS[mod1(series_idx, length(FIG7_COLORS))]
        offset = config["fig7_offset_x_axis"] * relative_x_offset(series_idx, n_series)

        for (pno_idx, pno) in enumerate(pn_orders)
            profile = profiles[series_id][pno]
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

function add_fig7_legend!(fig::Figure, target_slot, series_labels)
    elements = [
        PolyElement(
            color=(FIG7_COLORS[mod1(idx, length(FIG7_COLORS))], VIOLIN_ALPHA),
            strokecolor=(FIG7_COLORS[mod1(idx, length(FIG7_COLORS))], VIOLIN_EDGE_ALPHA),
            strokewidth=1.5,
        )
        for idx in eachindex(series_labels)
    ]

    Legend(target_slot, elements, series_labels;
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

function build_figure(config::Dict, profiles, series_ids, series_labels)
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
            panel_limits(panel_orders, profiles, series_ids) :
            config["fig7_y_limits"][panel_idx]

        colsize!(fig.layout, panel_idx, Relative(length(panel_orders) / total_width))
        build_top_label_row!(top_grid, panel_orders)
        configure_panel_axis!(ax, panel_orders, y_limits;
            ylabel=L"\delta\varphi_p",
            show_ylabel=panel_idx == 1)
        plot_panel!(ax, panel_orders, profiles, config, series_ids)
    end

    rowsize!(fig.layout, 1, Relative(FIG7_TOP_ROW_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG7_TOP_ROW_FRACTION))
    colgap!(fig.layout, FIG7_COL_GAP)
    rowgap!(fig.layout, FIG7_LABEL_ROW_GAP)
    add_fig7_legend!(fig, fig[2, n_panels], series_labels)

    return fig
end

function run_plot_fig7(config::Dict)
    sources = fig7_plot_sources(config)
    profiles = Dict{String, Dict{String, NamedTuple}}()
    series_ids = String[]
    series_labels = Any[]

    for (idx, source) in enumerate(sources)
        source_config = source.config
        use_realization = source_config["fig7_number_of_events_single_realization"] !== nothing

        series_id = "series_$(idx)"
        selected_dphi_k, selected_delta_k = if use_realization
            fisher_results_file = get_fisher_results_file(source_config)
            isfile(fisher_results_file) ||
                throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))
            read_realization_likelihood_data_from_fisher(
                fisher_results_file,
                source_config,
                source.waveform_family,
            )
        else
            population_results_file = get_population_results_file(source_config)
            isfile(population_results_file) ||
                throw(ArgumentError("Population results file does not exist: $(population_results_file)"))
            read_selected_likelihood_data(
                population_results_file,
                source_config,
                source.waveform_family,
            )
        end
        profiles[series_id] = build_posterior_profiles(
            selected_dphi_k,
            selected_delta_k,
            source_config,
            source.label,
        )
        push!(series_ids, series_id)
        push!(series_labels, source.label)
    end

    fig = build_figure(config, profiles, series_ids, series_labels)

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

### leave unchanged
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))
####

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
