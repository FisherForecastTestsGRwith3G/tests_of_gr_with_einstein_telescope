using TOML
using HDF5
using CairoMakie
using LaTeXStrings
using Colors
using Printf

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

const FIG8_PN_ORDERS = ["-1", "0", "0.5", "1", "1.5", "2", "3", "3.5", "log(2.5)", "log(3.)"]
const FIG8_SIZE = (1500, 560)
const FIG8_FIGURE_PADDING = (18, 36, 14, 10)
const FIG8_AXIS_LABELSIZE = 30
const FIG8_TICK_LABELSIZE = 22
const FIG8_LEGEND_LABELSIZE = 22
const FIG8_LINEWIDTH = 3.0
const FIG8_REFERENCE_LINEWIDTH = 2.5
const FIG8_CONNECTOR_LINEWIDTH = 1.4
const FIG8_MINUS_ONE_PN_SCALE = 1000.0
const FIG8_INSET_AXIS_SCALE = 10.0
const FIG8_CI = 0.90

const FIG8_MAIN_LIMITS = (
    x=(-0.010, 0.020),
    y=(0.0, 0.050),
)
const FIG8_INSET_LIMITS = (
    x=(-0.001, 0.001),
    y=(0.0, 0.002),
)

const FIG8_PN_LABELS = Dict(
    "-1" => L"1000\times\varphi_{-2}",
    "0" => L"\varphi_{0}",
    "0.5" => L"\varphi_{1}",
    "1" => L"\varphi_{2}",
    "1.5" => L"\varphi_{3}",
    "2" => L"\varphi_{4}",
    "log(2.5)" => L"\varphi_{5\ell}",
    "3" => L"\varphi_{6}",
    "log(3.)" => L"\varphi_{6\ell}",
    "3.5" => L"\varphi_{7}",
)

const FIG8_PN_COLORS = Dict(
    "-1" => "#3B4CC0",
    "0" => "#5F7FE8",
    "0.5" => "#86A9FC",
    "1" => "#ADC9FD",
    "1.5" => "#CFDAEA",
    "2" => "#E8D3C5",
    "3" => "#F2A17F",
    "3.5" => "#E66C53",
    "log(2.5)" => "#D1493F",
    "log(3.)" => "#B40426",
)

const FIG8_DEFAULT_MU_LIMITS = Dict{String, Union{Nothing, Tuple{Float64, Float64}}}(
    "-1" => nothing,
    "0" => (-0.002, 0.002),
    "0.5" => (-0.002, 0.002),
    "1" => (-0.002, 0.002),
    "1.5" => (-0.002, 0.002),
    "2" => (-0.020, 0.020),
    "log(2.5)" => (-0.020, 0.020),
    "3" => (-0.020, 0.020),
    "log(3.)" => (-0.040, 0.040),
    "3.5" => (-0.020, 0.020),
)
const FIG8_DEFAULT_SIGMA_LIMITS = Dict{String, Union{Nothing, Tuple{Float64, Float64}}}(
    "-1" => nothing,
    "0" => (0.0, 0.003),
    "0.5" => (0.0, 0.003),
    "1" => (0.0, 0.003),
    "1.5" => (0.0, 0.003),
    "2" => (0.0, 0.050),
    "log(2.5)" => (0.0, 0.050),
    "3" => (0.0, 0.050),
    "log(3.)" => (0.0, 0.060),
    "3.5" => (0.0, 0.050),
)

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_fig8_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"], "fig_8")
    stem = "fig8_$(config["plot_tag"])"
    return (
        joinpath(output_dir, stem * ".png"),
        joinpath(output_dir, stem * ".pdf"),
    )
end

function read_limit_tuple(raw_limits)
    limits = Float64.(raw_limits)
    length(limits) == 2 || throw(ArgumentError("Figure 8 axis limits must have two values."))
    limits[1] < limits[2] || throw(ArgumentError("Figure 8 lower limit must be smaller than upper limit."))
    return (limits[1], limits[2])
end

function read_optional_pn_limits(raw_fig8_config::Dict, key::AbstractString, defaults)
    limits = copy(defaults)
    haskey(raw_fig8_config, key) || return limits

    raw_limits = raw_fig8_config[key]
    raw_limits isa AbstractDict ||
        throw(ArgumentError("`plots.fig8.$(key)` must be a table keyed by PN order."))

    for (pno, pno_limits) in raw_limits
        pno in keys(limits) || throw(ArgumentError("Unknown PN order `$(pno)` in `plots.fig8.$(key)`."))
        limits[String(pno)] = isempty(pno_limits) ? nothing : read_limit_tuple(pno_limits)
    end

    return limits
end

function read_fig8_config(config_file::AbstractString)
    config_file = abspath(config_file)
    config = read_config(config_file)
    raw_config = TOML.parsefile(config_file)
    fig8_config = get(get(raw_config, "plots", Dict{String, Any}()), "fig8", Dict{String, Any}())

    config["fig8_config_file"] = config_file
    config["fig8_waveform_family"] = String(get(fig8_config, "waveform_family", "PhenomHM"))
    config["fig8_n_events"] = Int(get(fig8_config, "n_events", config["n_events"]))
    config["fig8_number_of_events_single_realization"] = Int(get(fig8_config, "number_of_events_single_realization", 10000))
    config["fig8_observed_events_direct"] = Bool(get(fig8_config, "n_events_and_number_of_events_single_realization_refer_directly_to_observed_events", false))
    config["fig8_realization_index"] = Int(get(fig8_config, "realization_index", 1))
    config["fig8_use_inspiral_snr_threshold"] = Bool(get(fig8_config, "use_inspiral_snr_threshold", false))
    config["fig8_grid_points"] = Int(get(fig8_config, "grid_points", 700))
    config["fig8_ci"] = Float64(get(fig8_config, "credible_interval", FIG8_CI))
    config["fig8_minus_one_pn_scale"] = Float64(get(fig8_config, "minus_one_pn_scale", FIG8_MINUS_ONE_PN_SCALE))
    config["fig8_inset_axis_scale"] = Float64(get(fig8_config, "inset_axis_scale", FIG8_INSET_AXIS_SCALE))
    config["fig8_main_x_limits"] = haskey(fig8_config, "main_x_limits") ? read_limit_tuple(fig8_config["main_x_limits"]) : FIG8_MAIN_LIMITS.x
    config["fig8_main_y_limits"] = haskey(fig8_config, "main_y_limits") ? read_limit_tuple(fig8_config["main_y_limits"]) : FIG8_MAIN_LIMITS.y
    config["fig8_inset_x_limits"] = haskey(fig8_config, "inset_x_limits") ? read_limit_tuple(fig8_config["inset_x_limits"]) : FIG8_INSET_LIMITS.x
    config["fig8_inset_y_limits"] = haskey(fig8_config, "inset_y_limits") ? read_limit_tuple(fig8_config["inset_y_limits"]) : FIG8_INSET_LIMITS.y
    config["fig8_pn_orders"] = haskey(fig8_config, "pn_orders") ? String.(fig8_config["pn_orders"]) : FIG8_PN_ORDERS
    config["fig8_mu_limits"] = read_optional_pn_limits(fig8_config, "mu_limits", FIG8_DEFAULT_MU_LIMITS)
    config["fig8_sigma_limits"] = read_optional_pn_limits(fig8_config, "sigma_limits", FIG8_DEFAULT_SIGMA_LIMITS)

    all(pno -> pno in config["pn_orders"], config["fig8_pn_orders"]) ||
        throw(ArgumentError("Every PN order in `plots.fig8.pn_orders` must also be present in `deviations.pn_orders`."))
    config["fig8_waveform_family"] in config["waveform_families"] ||
        throw(ArgumentError("Waveform family `$(config["fig8_waveform_family"])` is not listed in `fisher.waveform`."))
    0.0 < config["fig8_ci"] < 1.0 ||
        throw(ArgumentError("`plots.fig8.credible_interval` must be between 0 and 1."))
    config["fig8_grid_points"] >= 10 ||
        throw(ArgumentError("`plots.fig8.grid_points` must be at least 10."))

    return config
end

function clean_likelihood_data(dphi_k::Vector{Float64}, delta_k::Vector{Float64})
    valid = isfinite.(dphi_k) .& isfinite.(delta_k) .& (delta_k .> 0.0)
    return dphi_k[valid], delta_k[valid]
end

function waveform_summary_selection_from_fisher(file, config::Dict, waveform_family::AbstractString)
    summary_selection = nothing
    for pno in config["fig8_pn_orders"]
        wf_group = file[createSED.pnoString(pno)][config["network"]][waveform_family]
        snr = Float64.(read(wf_group, "snr"))
        isnr = Float64.(read(wf_group, "isnr"))
        invc = Bool.(read(wf_group, "invc"))
        delta_k = Float64.(read(wf_group, "delta_k"))

        pno_selection = BitVector(
            (snr .> config["snr_threshold"]) .&
            invc .&
            isfinite.(delta_k) .&
            (delta_k .> 0.0)
        )
        if config["fig8_use_inspiral_snr_threshold"]
            pno_selection .&= (isnr .> config["snr_inspiral_threshold"])
        end
        summary_selection = isnothing(summary_selection) ? pno_selection : BitVector(summary_selection .& pno_selection)
    end
    return summary_selection
end

function realization_indices(summary_selection::BitVector, n_events_total::Integer, config::Dict)
    events_per_realization = config["fig8_number_of_events_single_realization"]
    n_events_requested = min(config["fig8_n_events"], n_events_total)
    observed_direct = config["fig8_observed_events_direct"]

    if observed_direct
        observed_indices = findall(summary_selection)
        n_events_used = min(n_events_requested, length(observed_indices))
        n_realizations = fld(n_events_used, events_per_realization)
        n_realizations > 0 ||
            throw(ArgumentError("Not enough observed events to build one Figure 8 realization."))
        realization_index = config["fig8_realization_index"]
        1 <= realization_index <= n_realizations ||
            throw(ArgumentError("Requested realization $(realization_index), but only $(n_realizations) are available."))
        first_idx = (realization_index - 1) * events_per_realization + 1
        last_idx = realization_index * events_per_realization
        chosen = observed_indices[first_idx:last_idx]
        println("Using observed-events realization $(realization_index)/$(n_realizations) with $(length(chosen)) observed events.")
        return chosen
    end

    n_realizations = fld(n_events_requested, events_per_realization)
    n_realizations > 0 ||
        throw(ArgumentError("Not enough catalog events to build one Figure 8 realization."))
    realization_index = config["fig8_realization_index"]
    1 <= realization_index <= n_realizations ||
        throw(ArgumentError("Requested realization $(realization_index), but only $(n_realizations) are available."))

    first_idx = (realization_index - 1) * events_per_realization + 1
    last_idx = realization_index * events_per_realization
    chunk_indices = collect(first_idx:last_idx)
    observed_chunk_indices = chunk_indices[summary_selection[chunk_indices]]
    println("Using catalog realization $(realization_index)/$(n_realizations): $(length(observed_chunk_indices)) observed events out of $(events_per_realization) catalog events.")
    return observed_chunk_indices
end

function read_fig8_likelihood_data(fisher_results_file::AbstractString, config::Dict)
    selected_dphi_k = Dict{String, Vector{Float64}}()
    selected_delta_k = Dict{String, Vector{Float64}}()
    waveform_family = config["fig8_waveform_family"]

    h5open(fisher_results_file, "r") do file
        reference_pno = first(config["fig8_pn_orders"])
        reference_group = file[createSED.pnoString(reference_pno)][config["network"]][waveform_family]
        n_events_total = length(read(reference_group, "dphi_k"))

        summary_selection = waveform_summary_selection_from_fisher(file, config, waveform_family)
        chosen_indices = realization_indices(summary_selection, n_events_total, config)
        isempty(chosen_indices) &&
            throw(ArgumentError("Figure 8 realization has no events after applying the observation cuts."))

        for pno in config["fig8_pn_orders"]
            wf_group = file[createSED.pnoString(pno)][config["network"]][waveform_family]
            dphi, delta = clean_likelihood_data(
                Float64.(read(wf_group, "dphi_k"))[chosen_indices],
                Float64.(read(wf_group, "delta_k"))[chosen_indices],
            )
            isempty(dphi) &&
                throw(ArgumentError("No finite selected likelihood data for waveform $(waveform_family), PN order $(pno)."))
            selected_dphi_k[pno] = dphi
            selected_delta_k[pno] = delta
        end
    end

    return selected_dphi_k, selected_delta_k
end

function automatic_grid_limits(dphi_k::Vector{Float64}, delta_k::Vector{Float64})
    hyperparam_dist = HierDist.hyperparamDistTIGER(dphi_k, delta_k)
    return HierDist.findOptimalGrid(hyperparam_dist; verbose=false)
end

function grid_limits_for_pno(config::Dict, pno::AbstractString, dphi_k::Vector{Float64}, delta_k::Vector{Float64})
    mu_limits = config["fig8_mu_limits"][pno]
    sigma_limits = config["fig8_sigma_limits"][pno]
    if isnothing(mu_limits) || isnothing(sigma_limits)
        auto_mu_min, auto_mu_max, auto_sigma_min, auto_sigma_max = automatic_grid_limits(dphi_k, delta_k)
        mu_limits = isnothing(mu_limits) ? (auto_mu_min, auto_mu_max) : mu_limits
        sigma_limits = isnothing(sigma_limits) ? (auto_sigma_min, auto_sigma_max) : sigma_limits
    end
    return mu_limits, sigma_limits
end

function build_fig8_contours(selected_dphi_k, selected_delta_k, config::Dict)
    contours = Dict{String, Vector{NamedTuple{(:mu, :sigma), Tuple{Vector{Float64}, Vector{Float64}}}}}()

    for pno in config["fig8_pn_orders"]
        println("Building $(round(Int, 100 * config["fig8_ci"]))% contour for PN order $(pno)")
        dphi_k = selected_dphi_k[pno]
        delta_k = selected_delta_k[pno]
        mu_limits, sigma_limits = grid_limits_for_pno(config, pno, dphi_k, delta_k)
        mu_grid = collect(range(mu_limits[1], mu_limits[2]; length=config["fig8_grid_points"]))
        sigma_grid = collect(range(max(0.0, sigma_limits[1]), sigma_limits[2]; length=config["fig8_grid_points"]))

        hyperparam_dist = HierDist.hyperparamDistTIGER(dphi_k, delta_k)
        p_mu_sigma, _, _, _, _ = HierDist.getDistributionOnGrid(mu_grid, sigma_grid, hyperparam_dist)
        mu_lines, sigma_lines = HierDist.getCredibleContourOGD(
            mu_grid,
            sigma_grid,
            p_mu_sigma,
            config["fig8_ci"],
        )
        contours[pno] = [
            (mu=mu_lines[idx], sigma=sigma_lines[idx])
            for idx in eachindex(mu_lines)
        ]
    end

    return contours
end

function fig8_line_style(pno::AbstractString)
    return pno == "-1" ? :dash : :solid
end

function fig8_scaled_contour(contour, pno::AbstractString, config::Dict)
    scale = pno == "-1" ? config["fig8_minus_one_pn_scale"] : 1.0
    return contour.mu .* scale, contour.sigma .* scale
end

function add_fig8_contours!(ax::Axis, contours, config::Dict; add_labels::Bool=false)
    for pno in config["fig8_pn_orders"]
        color = parse(Colorant, FIG8_PN_COLORS[pno])
        for contour in contours[pno]
            mu_values, sigma_values = fig8_scaled_contour(contour, pno, config)
            lines!(
                ax,
                mu_values,
                sigma_values;
                color=color,
                linewidth=FIG8_LINEWIDTH,
                linestyle=fig8_line_style(pno),
                label=add_labels ? FIG8_PN_LABELS[pno] : nothing,
            )
        end
    end
    return ax
end

function scaled_ticklabels(ticks::Vector{Float64}, scale::Float64; digits::Int=3, blank_indices=Int[])
    labels = String[]
    for (idx, tick) in enumerate(ticks)
        if idx in blank_indices
            push!(labels, "")
        else
            value = tick * scale
            abs(value) < 0.5 * 10.0^(-digits) && (value = 0.0)
            push!(labels, @sprintf("%.*f", digits, value))
        end
    end
    return labels
end

function configure_main_axis!(ax::Axis, config::Dict)
    ax.xlabel = L"\mu"
    ax.ylabel = L"\sigma"
    ax.xlabelsize = FIG8_AXIS_LABELSIZE
    ax.ylabelsize = FIG8_AXIS_LABELSIZE
    ax.xticklabelsize = FIG8_TICK_LABELSIZE
    ax.yticklabelsize = FIG8_TICK_LABELSIZE
    ax.xtickalign = 1
    ax.ytickalign = 1
    ax.xgridvisible = false
    ax.ygridvisible = false
    ax.xticks = -0.010:0.005:0.020
    ax.yticks = 0.0:0.01:0.05
    xlims!(ax, config["fig8_main_x_limits"]...)
    ylims!(ax, config["fig8_main_y_limits"]...)
    vlines!(ax, [0.0]; color=(:gray45, 0.45), linestyle=:dash, linewidth=FIG8_REFERENCE_LINEWIDTH)
    return ax
end

function configure_inset_axis!(ax::Axis, config::Dict)
    inset_scale = config["fig8_inset_axis_scale"]
    ax.xlabel = latexstring("$(round(Int, inset_scale))\\times\\mu")
    ax.ylabel = latexstring("$(round(Int, inset_scale))\\times\\sigma")
    ax.yaxisposition = :right
    ax.xlabelsize = FIG8_AXIS_LABELSIZE
    ax.ylabelsize = FIG8_AXIS_LABELSIZE
    ax.xticklabelsize = FIG8_TICK_LABELSIZE
    ax.yticklabelsize = FIG8_TICK_LABELSIZE
    ax.xtickalign = 1
    ax.ytickalign = 1
    ax.xgridvisible = false
    ax.ygridvisible = false
    xlims!(ax, config["fig8_inset_x_limits"]...)
    ylims!(ax, config["fig8_inset_y_limits"]...)

    x_ticks = collect(range(config["fig8_inset_x_limits"][1], config["fig8_inset_x_limits"][2]; length=9))
    y_ticks = collect(range(config["fig8_inset_y_limits"][1], config["fig8_inset_y_limits"][2]; length=5))
    ax.xticks = (x_ticks, scaled_ticklabels(x_ticks, inset_scale; digits=3, blank_indices=[1, 2, 4, 6, 8]))
    ax.yticks = (y_ticks, scaled_ticklabels(y_ticks, inset_scale; digits=2, blank_indices=[2, 4]))

    vlines!(ax, [0.0]; color=(:gray45, 0.45), linestyle=:dash, linewidth=FIG8_REFERENCE_LINEWIDTH)
    return ax
end

function add_zoom_box!(ax::Axis, config::Dict)
    x1, x2 = config["fig8_inset_x_limits"]
    y1, y2 = config["fig8_inset_y_limits"]
    rect_x = [x1, x2, x2, x1, x1]
    rect_y = [y1, y1, y2, y2, y1]
    lines!(ax, rect_x, rect_y; color=:black, linewidth=FIG8_CONNECTOR_LINEWIDTH)
    return ax
end

function axis_data_to_figure_pixel(ax::Axis, x::Real, y::Real)
    local_point = Makie.project(ax.scene, :data, :pixel, Point2f(x, y))
    viewport = ax.scene.viewport[]
    return Point2f(
        viewport.origin[1] + local_point[1],
        viewport.origin[2] + local_point[2],
    )
end

function add_zoom_connectors!(fig::Figure, ax::Axis, ax_inset::Axis, config::Dict)
    resize_to_layout!(fig)

    x1, x2 = config["fig8_inset_x_limits"]
    _, y2 = config["fig8_inset_y_limits"]
    zoom_top_left = axis_data_to_figure_pixel(ax, x1, y2)
    zoom_top_right = axis_data_to_figure_pixel(ax, x2, y2)

    inset_viewport = ax_inset.scene.viewport[]
    inset_left = Float32(inset_viewport.origin[1])
    inset_bottom = Float32(inset_viewport.origin[2])
    inset_top = Float32(inset_viewport.origin[2] + inset_viewport.widths[2])

    lines!(
        fig.scene,
        [zoom_top_left[1], inset_left],
        [zoom_top_left[2], inset_top];
        color=:black,
        linestyle=:dash,
        linewidth=FIG8_CONNECTOR_LINEWIDTH,
        space=:pixel,
    )
    lines!(
        fig.scene,
        [zoom_top_right[1], inset_left],
        [zoom_top_right[2], inset_bottom];
        color=:black,
        linestyle=:dash,
        linewidth=FIG8_CONNECTOR_LINEWIDTH,
        space=:pixel,
    )
    return fig
end

function add_fig8_legend!(fig::Figure, target_slot, config::Dict)
    elements = [
        LineElement(
            color=parse(Colorant, FIG8_PN_COLORS[pno]),
            linewidth=FIG8_LINEWIDTH,
            linestyle=fig8_line_style(pno),
        )
        for pno in config["fig8_pn_orders"]
    ]
    labels = [FIG8_PN_LABELS[pno] for pno in config["fig8_pn_orders"]]

    Legend(
        target_slot,
        elements,
        labels;
        framevisible=false,
        orientation=:horizontal,
        nbanks=2,
        tellwidth=false,
        tellheight=true,
        labelsize=FIG8_LEGEND_LABELSIZE,
        patchsize=(48, 18),
        colgap=24,
        rowgap=4,
    )
    return fig
end

function build_fig8_figure(contours, config::Dict)
    CairoMakie.activate!()
    fig = Figure(size=FIG8_SIZE, backgroundcolor=:white, figure_padding=FIG8_FIGURE_PADDING)
    legend_grid = GridLayout(fig[1, 1:2])
    ax = Axis(fig[2, 1], backgroundcolor=:white)
    inset_grid = GridLayout(fig[2, 2])
    ax_inset = Axis(inset_grid[1, 1], backgroundcolor=:white)
    Box(inset_grid[2, 1]; color=(:white, 0.0), strokecolor=(:white, 0.0))

    colsize!(fig.layout, 1, Relative(0.56))
    colsize!(fig.layout, 2, Relative(0.44))
    rowsize!(fig.layout, 1, Auto(0.16))
    rowsize!(fig.layout, 2, Relative(0.84))
    colgap!(fig.layout, -48)
    rowgap!(fig.layout, 0)
    rowsize!(inset_grid, 1, Relative(0.78))
    rowsize!(inset_grid, 2, Relative(0.22))
    rowgap!(inset_grid, 0)

    configure_main_axis!(ax, config)
    configure_inset_axis!(ax_inset, config)
    add_fig8_contours!(ax, contours, config)
    add_fig8_contours!(ax_inset, contours, config)
    add_zoom_box!(ax, config)
    add_fig8_legend!(fig, legend_grid[1, 1], config)
    add_zoom_connectors!(fig, ax, ax_inset, config)

    return fig
end

function run_plot_fig8(config::Dict)
    fisher_results_file = get_fisher_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    println("Using Fisher results file: $(fisher_results_file)")
    println("Waveform family: $(config["fig8_waveform_family"])")
    selected_dphi_k, selected_delta_k = read_fig8_likelihood_data(fisher_results_file, config)
    contours = build_fig8_contours(selected_dphi_k, selected_delta_k, config)
    fig = build_fig8_figure(contours, config)

    png_output_file, pdf_output_file = get_fig8_output_files(config)
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
    config = read_fig8_config(config_file)
    return run_plot_fig8(config)
end

### leave unchanged
fontsize_theme = Theme(fontsize=24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts=(
    regular=joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold=joinpath(mt_fonts_dir, "NewCM10-Bold.otf"),
))
####

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
