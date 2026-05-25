using TOML
using HDF5
using CairoMakie
using LaTeXStrings
using Statistics
import JSON

## WARNING: NOT HUMAN CONTROLLED YET! 

include("_config_parser.jl")
include("_plot_style.jl")
include("../create_single_event_datasets/createSED.jl")

const PHENOM_HM = "PhenomHM"
const DETECTOR_ORDER = ["ET_0_15km", "ET_45_15km", "ETS"]
const DETECTOR_LABELS = Dict(
    "ET_0_15km" => L"\mathrm{ET}\ \mathrm{2L\_0}",
    "ET_45_15km" => L"\mathrm{ET}\ \mathrm{2L\_45}",
    "ETS" => L"\mathrm{ET}\ \Delta",
)
const DETECTOR_COLORS = Dict(
    "ET_0_15km" => FIG9_IMPROVEMENT_HIGH_COLOR,
    "ET_45_15km" => FIG9_IMPROVEMENT_MIDDLE_COLOR,
    "ETS" => FIG9_IMPROVEMENT_LOW_COLOR,
)
const DETECTOR_FILL_COLORS = Dict(
    "ET_0_15km" => FIG9_LIGHT_IMPROVEMENT_HIGH_COLOR,
    "ET_45_15km" => FIG9_LIGHT_IMPROVEMENT_MIDDLE_COLOR,
    "ETS" => FIG9_LIGHT_IMPROVEMENT_LOW_COLOR,
)
const GWTC3_COLOR = :black
const GUIDE_FONT_SIZE = 32
const TICK_FONT_SIZE = 32
const LEGEND_FONT_SIZE = 28
const FIG2_SIZE = (1800, 900)
const FIG2_TOP_ROW_FRACTION = 0.11
const FIG2_COL_GAP = 20
const FIG2_LABEL_ROW_GAP = 10
const HIST_MAX_WIDTH = 0.20
const SUMMARY_X_OFFSET = -0.0455
const RIGHT_PANEL_Y_LIMITS = (5e-6, 5e2)
const DETECTOR_OFFSETS = Dict(
    "ET_0_15km" => -0.13,
    "ET_45_15km" => 0.0,
    "ETS" => 0.13,
)
const GWTC3_REFERENCE_FILE = joinpath(@__DIR__, "lvk_gwtc_3_results_2025.json")

function load_gwtc3_reference()
    all_results = JSON.parsefile(GWTC3_REFERENCE_FILE; dicttype=Dict{String, Any})
    reference = get(all_results, "GWTC-3 (SEOB)", nothing)
    reference === nothing && throw(ArgumentError("Missing `GWTC-3 (SEOB)` in $(GWTC3_REFERENCE_FILE)"))
    return Dict{String, Float64}(key => Float64(value) for (key, value) in reference)
end

const GWTC3_REFERENCE = load_gwtc3_reference()

function get_population_results_file(config::Dict)
    return joinpath(
        @__DIR__,
        config["bootstrap_outdir"],
        "population_results_$(config["bootstrap_tag"]).h5"
    )
end

function get_plot_output_files(configs::Vector{Dict{String, Any}})
    first_config = configs[1]
    output_dir = joinpath(@__DIR__, first_config["plot_outdir"], "fig_2")
    plot_tag = "detector_comparison_$(first_config["plot_tag"])"
    return (
        joinpath(output_dir, "fig2_$(plot_tag).png"),
        joinpath(output_dir, "fig2_$(plot_tag).pdf"),
    )
end

function get_histogram_edges(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)
    positive_samples = samples[samples .> 0.0]
    isempty(positive_samples) && return nothing

    y_min = minimum(positive_samples)
    y_max = maximum(positive_samples)
    if y_min == y_max
        return exp10.(range(log10(y_min) - 0.25, log10(y_max) + 0.25, length=11))
    end

    if isnothing(n_bins)
        n_bins = clamp(round(Int, sqrt(length(positive_samples))), 8, 40)
    end

    return exp10.(range(log10(y_min), log10(y_max), length=n_bins + 1))
end

function histogram_profile(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)
    edges = get_histogram_edges(samples; n_bins=n_bins)
    edges === nothing && return nothing

    counts = zeros(Int, length(edges) - 1)
    for value in samples
        value <= 0.0 && continue
        idx = searchsortedlast(edges, value)
        idx = clamp(idx, 1, length(edges) - 1)
        value == edges[end] && (idx = length(edges) - 1)
        counts[idx] += 1
    end

    return edges, counts
end

function mirrored_histogram_bars!(ax::Axis, x0::Real, samples::Vector{Float64}, fill_color, edge_color;
    max_width::Float64=HIST_MAX_WIDTH, n_bins::Union{Nothing, Int}=nothing)
    profile = histogram_profile(samples; n_bins=n_bins)
    profile === nothing && return nothing

    edges, counts = profile
    max_count = maximum(counts)
    max_count == 0 && return nothing

    widths = max_width .* (counts ./ max_count)
    for idx in eachindex(counts)
        widths[idx] == 0 && continue
        y1 = edges[idx]
        y2 = edges[idx + 1]
        poly!(ax, Point2f[
                (x0 - widths[idx], y1),
                (x0 + widths[idx], y1),
                (x0 + widths[idx], y2),
                (x0 - widths[idx], y2),
            ];
            color=fill_color, strokecolor=:transparent)
    end

    left_x = Float64[]
    left_y = Float64[]
    right_x = Float64[]
    right_y = Float64[]
    for idx in eachindex(counts)
        push!(left_y, edges[idx], edges[idx + 1])
        push!(right_y, edges[idx], edges[idx + 1])
        push!(left_x, x0 - widths[idx], x0 - widths[idx])
        push!(right_x, x0 + widths[idx], x0 + widths[idx])
    end
    lines!(ax, left_x, left_y; color=edge_color, linewidth=2)
    lines!(ax, right_x, right_y; color=edge_color, linewidth=2)

    return nothing
end

function add_bootstrap_summary!(ax::Axis, x0::Real, samples::Vector{Float64}, color)
    positive_samples = samples[samples .> 0.0]
    isempty(positive_samples) && return nothing

    q05 = quantile(positive_samples, 0.05)
    q50 = quantile(positive_samples, 0.50)
    q95 = quantile(positive_samples, 0.95)
    x_marker = x0 + SUMMARY_X_OFFSET

    lines!(ax, [x_marker, x_marker], [q05, q95]; color=color, linewidth=3)
    scatter!(ax, [x_marker], [q50]; color=color, markersize=15, strokecolor=:white, strokewidth=1.0)
    return nothing
end

function add_gwtc3_reference!(ax::Axis, x0::Real, pno::String)
    haskey(GWTC3_REFERENCE, pno) || return nothing
    scatter!(ax, [x0], [GWTC3_REFERENCE[pno]];
        color=GWTC3_COLOR, marker=:diamond, markersize=18,
        strokecolor=:white, strokewidth=1.0)
    return nothing
end

function read_detector_plot_data(population_results_file::AbstractString, config::Dict)
    single_event_constraints = Dict{String, Vector{Float64}}()
    bootstrap_constraints = Dict{String, Vector{Float64}}()

    h5open(population_results_file, "r") do file
        network_group = file[config["network"]]
        haskey(network_group, PHENOM_HM) ||
            throw(ArgumentError("Missing waveform group `$(PHENOM_HM)` in $(population_results_file) for detector $(config["network"])."))
        waveform_group = network_group[PHENOM_HM]

        for pno in config["pn_orders"]
            pno_group = waveform_group[createSED.pnoString(pno)]
            single_event_constraints[pno] = read(pno_group, "single_event_constraints")
            bootstrap_constraints[pno] = read(pno_group, "bootstrap_constraints")
        end
    end

    return single_event_constraints, bootstrap_constraints
end

function read_all_plot_data(configs::Vector{Dict{String, Any}})
    single_event_by_detector = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_by_detector = Dict{String, Dict{String, Vector{Float64}}}()

    for config in configs
        population_results_file = get_population_results_file(config)
        isfile(population_results_file) ||
            throw(ArgumentError("Population results file does not exist: $(population_results_file)"))

        single_event, bootstrap = read_detector_plot_data(population_results_file, config)
        single_event_by_detector[config["network"]] = single_event
        bootstrap_by_detector[config["network"]] = bootstrap
    end

    return single_event_by_detector, bootstrap_by_detector
end

function histogram_limits(
    pn_orders::AbstractVector{<:AbstractString},
    single_event_by_detector::Dict{String, Dict{String, Vector{Float64}}},
)
    y_min = Inf
    y_max = 0.0

    for pno in pn_orders
        for detector in DETECTOR_ORDER
            samples = get(get(single_event_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            edges = get_histogram_edges(samples)
            edges === nothing && continue
            y_min = min(y_min, first(edges))
            y_max = max(y_max, last(edges))
        end
    end

    isfinite(y_min) || return nothing
    return (y_min, y_max)
end

function panel_data_limits(
    pn_orders::AbstractVector{<:AbstractString},
    single_event_by_detector::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_by_detector::Dict{String, Dict{String, Vector{Float64}}},
)
    values = Float64[]

    for pno in pn_orders
        for detector in DETECTOR_ORDER
            single_event = get(get(single_event_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            append!(values, single_event[single_event .> 0.0])

            bootstrap = get(get(bootstrap_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            positive_bootstrap = bootstrap[bootstrap .> 0.0]
            isempty(positive_bootstrap) || append!(values, quantile(positive_bootstrap, [0.05, 0.50, 0.95]))
        end

        if haskey(GWTC3_REFERENCE, pno) && GWTC3_REFERENCE[pno] > 0.0
            push!(values, GWTC3_REFERENCE[pno])
        end
    end

    isempty(values) && return nothing
    return (minimum(values), maximum(values))
end

function expand_log_limits(y_limits::Tuple{<:Real, <:Real}; lower_pad_decades::Real=0.5, upper_pad_decades::Real=0.08)
    y_min, y_max = y_limits
    log_y_min = log10(y_min)
    log_y_max = log10(y_max)
    return (
        exp10(log_y_min - lower_pad_decades),
        exp10(log_y_max + upper_pad_decades),
    )
end

function decade_ticks(y_limits::Tuple{<:Real, <:Real})
    y_min, y_max = y_limits
    decade_min = floor(Int, log10(y_min))
    decade_max = ceil(Int, log10(y_max))
    exponents = collect(decade_min:decade_max)
    length(exponents) > 2 && (exponents = exponents[2:end-1])
    values = exp10.(exponents)
    labels = [latexstring("10^{", exponent, "}") for exponent in exponents]
    return values, labels
end

function top_pno_label(pno::AbstractString)
    if startswith(pno, "log(") && endswith(pno, ")")
        order = chop(chop(pno; head=4); tail=1)
        endswith(order, ".") && (order = chop(order; tail=1))
        return latexstring(order, raw"\,\mathrm{PN}^{(\ell)}")
    end

    return latexstring(pno, raw"\,\mathrm{PN}")
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

function configure_panel_axis!(ax::Axis, pn_orders::AbstractVector{<:AbstractString}, y_limits::Tuple{<:Real, <:Real}, y_tick_spec;
    ylabel="", show_ylabel::Bool=true)
    ax.xlabel = "PN order"
    ax.ylabel = show_ylabel ? ylabel : ""
    ax.xlabelsize = GUIDE_FONT_SIZE
    ax.ylabelsize = GUIDE_FONT_SIZE
    ax.xticks = (collect(1:length(pn_orders)), createSED.pnoLatex.(pn_orders))
    ax.yticks = y_tick_spec
    ax.xticklabelsize = TICK_FONT_SIZE
    ax.yticklabelsize = TICK_FONT_SIZE
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = false
    ax.yminorgridvisible = true
    ax.yminorgridcolor = (:gray70, 0.35)
    ax.yminorticks = IntervalsBetween(9)
    ax.xtickalign = 1
    ax.ytickalign = 1
    xlims!(ax, 0.5, length(pn_orders) + 0.5)
    ylims!(ax, y_limits...)
    return ax
end

function build_panel!(ax::Axis, pn_orders::AbstractVector{<:AbstractString},
    single_event_by_detector::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_by_detector::Dict{String, Dict{String, Vector{Float64}}})
    for (idx, pno) in enumerate(pn_orders)
        for detector in DETECTOR_ORDER
            x0 = idx + DETECTOR_OFFSETS[detector]
            color = DETECTOR_COLORS[detector]
            fill_color = DETECTOR_FILL_COLORS[detector]
            data = get(get(single_event_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            boot = get(get(bootstrap_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            mirrored_histogram_bars!(ax, x0, data, fill_color, color)
            add_bootstrap_summary!(ax, x0, boot, color)
        end
        add_gwtc3_reference!(ax, idx, pno)
    end
    return ax
end

function add_fig2_legend!(fig::Figure, target_slot)
    summary_elements = [
        MarkerElement(color=:black, marker=:circle, markersize=15, strokecolor=:white, strokewidth=1.0),
        LineElement(color=:black, linewidth=3),
        MarkerElement(color=GWTC3_COLOR, marker=:diamond, markersize=18, strokecolor=:white, strokewidth=1.0),
    ]
    summary_labels = [
        L"\text{Population constraint}",
        L"\text{Population constraint 90\% CI}",
        L"\text{GWTC-3 TGR (SEOBNRv4\_ROM)}",
    ]

    elements = Any[]
    labels = Any[]
    for idx in eachindex(DETECTOR_ORDER)
        detector = DETECTOR_ORDER[idx]
        push!(elements, PolyElement(color=DETECTOR_FILL_COLORS[detector], strokecolor=DETECTOR_COLORS[detector]))
        push!(labels, DETECTOR_LABELS[detector])
        push!(elements, summary_elements[idx])
        push!(labels, summary_labels[idx])
    end

    Legend(target_slot, elements, labels;
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:bottom,
        margin=(10, 10, 10, 10),
        framevisible=true,
        backgroundcolor=(:white, 0.9),
        nbanks=2,
        colgap=24,
        rowgap=8,
        labelsize=LEGEND_FONT_SIZE)
    return fig
end

function build_figure(
    left_orders::Vector{String},
    right_orders::Vector{String},
    single_event_by_detector::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_by_detector::Dict{String, Dict{String, Vector{Float64}}},
    left_y_limits::Tuple{<:Real, <:Real},
    right_y_limits::Tuple{<:Real, <:Real},
)
    CairoMakie.activate!()

    fig = Figure(size=FIG2_SIZE, backgroundcolor=:white)
    left_top_grid = GridLayout(fig[1, 1])
    right_top_grid = GridLayout(fig[1, 2])
    left_ax = Axis(fig[2, 1], backgroundcolor=:white, yscale=log10)
    right_ax = Axis(fig[2, 2], backgroundcolor=:white, yscale=log10)

    panel_width_weights = Float64[length(left_orders), length(right_orders)]
    width_total = sum(panel_width_weights)
    colsize!(fig.layout, 1, Relative(panel_width_weights[1] / width_total))
    colsize!(fig.layout, 2, Relative(panel_width_weights[2] / width_total))
    rowsize!(fig.layout, 1, Relative(FIG2_TOP_ROW_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG2_TOP_ROW_FRACTION))
    colgap!(fig.layout, FIG2_COL_GAP)
    rowgap!(fig.layout, FIG2_LABEL_ROW_GAP)

    build_top_label_row!(left_top_grid, left_orders)
    build_top_label_row!(right_top_grid, right_orders)

    configure_panel_axis!(left_ax, left_orders, left_y_limits, decade_ticks(left_y_limits);
        ylabel=L"|\delta\varphi_p|", show_ylabel=true)
    configure_panel_axis!(right_ax, right_orders, right_y_limits, decade_ticks(right_y_limits);
        ylabel="", show_ylabel=false)

    build_panel!(left_ax, left_orders, single_event_by_detector, bootstrap_by_detector)
    build_panel!(right_ax, right_orders, single_event_by_detector, bootstrap_by_detector)
    add_fig2_legend!(fig, fig[2, 2])

    return fig
end

function assert_same_vector(configs::Vector{Dict{String, Any}}, key::String)
    reference = configs[1][key]
    for config in configs[2:end]
        config[key] == reference ||
            throw(ArgumentError("Config mismatch: `$(key)` must agree across the three configs."))
    end
end

function validate_detector_configs(configs::Vector{Dict{String, Any}})
    length(configs) == length(DETECTOR_ORDER) ||
        throw(ArgumentError("Expected exactly three configs: ET_0_15km, ET_45_15km, ETS."))

    for (idx, expected_detector) in enumerate(DETECTOR_ORDER)
        actual_detector = configs[idx]["network"]
        actual_detector == expected_detector ||
            throw(ArgumentError("Config $(idx) must use detector `$(expected_detector)`, got `$(actual_detector)`."))

        PHENOM_HM in configs[idx]["waveform_families"] ||
            throw(ArgumentError("Config $(idx) must include `$(PHENOM_HM)` in `fisher.waveform`."))
    end

    for key in ("n_catalog", "n_sample", "pn_orders", "mu", "sigma")
        assert_same_vector(configs, key)
    end

    for key in ("snr_threshold", "snr_inspiral_threshold", "select_before_bootstrap")
        assert_same_vector(configs, key)
    end

    return configs
end

function run_plot_fig2(configs::Vector{Dict{String, Any}})
    validate_detector_configs(configs)
    single_event_by_detector, bootstrap_by_detector = read_all_plot_data(configs)

    pn_orders = configs[1]["pn_orders"]
    left_orders = ["-1"]
    right_orders = filter(pno -> pno != "-1", pn_orders)

    bootstrap_medians = Dict{String, Float64}()
    for pno in ("-1", "0")
        positive_samples = Float64[]
        for detector in DETECTOR_ORDER
            samples = get(get(bootstrap_by_detector, detector, Dict{String, Vector{Float64}}()), pno, Float64[])
            append!(positive_samples, samples[samples .> 0.0])
        end
        isempty(positive_samples) || (bootstrap_medians[pno] = median(positive_samples))
    end

    right_scale = 1.0
    if haskey(bootstrap_medians, "-1") && haskey(bootstrap_medians, "0")
        median_ratio = bootstrap_medians["-1"] / bootstrap_medians["0"]
        if isfinite(median_ratio) && median_ratio > 0.0
            right_scale = 10.0 ^ round(Int, log10(median_ratio))
        end
    end

    right_y_limits = RIGHT_PANEL_Y_LIMITS
    left_y_limits = (right_y_limits[1] * right_scale, right_y_limits[2] * right_scale)

    fig = build_figure(
        left_orders,
        right_orders,
        single_event_by_detector,
        bootstrap_by_detector,
        left_y_limits,
        right_y_limits,
    )

    png_output_file, pdf_output_file = get_plot_output_files(configs)
    mkpath(dirname(png_output_file))
    save(png_output_file, fig)
    save(pdf_output_file, fig)
    println("Saved figure to $(png_output_file)")
    println("Saved figure to $(pdf_output_file)")

    return png_output_file
end

function main(args=ARGS)
    if length(args) != 3
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <ET_0_15km_config.toml> <ET_45_15km_config.toml> <ETS_config.toml>"))
    end

    config_files = abspath.(args)
    for config_file in config_files
        println("Using config file: $(config_file)")
    end

    configs = Dict{String, Any}[read_config(config_file) for config_file in config_files]
    return run_plot_fig2(configs)
end

fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
