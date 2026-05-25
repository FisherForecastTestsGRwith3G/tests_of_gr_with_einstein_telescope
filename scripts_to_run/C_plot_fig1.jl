using TOML
using HDF5
using CairoMakie
using LaTeXStrings
using Statistics
import JSON

include("_config_parser.jl")
include("_plot_style.jl")
include("../create_single_event_datasets/createSED.jl")

const PHENOM_D_COLOR  = FIG9_IMPROVEMENT_LOW_COLOR
const PHENOM_HM_COLOR = FIG9_IMPROVEMENT_HIGH_COLOR
const PHENOM_D_FILL_COLOR  = FIG9_LIGHT_IMPROVEMENT_LOW_COLOR
const PHENOM_HM_FILL_COLOR = FIG9_LIGHT_IMPROVEMENT_HIGH_COLOR
const GWTC3_COLOR     = :black
const TITLE_FONT_SIZE = 28
const GUIDE_FONT_SIZE = 32
const TICK_FONT_SIZE  = 32
const LEGEND_FONT_SIZE = 28
const FIG1_SIZE = (1800, 760)
const FIG1_TOP_ROW_FRACTION = 0.11
const FIG1_COL_GAP = 20
const FIG1_LABEL_ROW_GAP = 10
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

function get_plot_output_file(config::Dict)
    return joinpath(
        @__DIR__,
        config["plot_outdir"],
        "fig1_$(config["plot_tag"]).pdf"
    )
end

function get_plot_output_files(config::Dict)
    configured_output_file = get_plot_output_file(config)
    output_dir = joinpath(dirname(configured_output_file), "fig_1")
    output_stem = splitext(basename(configured_output_file))[1]
    return (
        joinpath(output_dir, output_stem * ".png"),
        joinpath(output_dir, output_stem * ".pdf"),
    )
end

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#
"""
    get_histogram_edges(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)

Return logarithmically spaced histogram bin edges for the positive entries in
`samples`. Returns `nothing` when no positive samples are present.
Lower (upper) bound of the bins is given by the smallest (largest) sample handed.
"""
function get_histogram_edges(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)
    positive_samples = samples[samples .> 0.0]
    if isempty(positive_samples)
        return nothing
    end

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

"""
    histogram_profile(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)

Build histogram bin counts for the positive entries in `samples` using logarithmically
spaced bin edges from [`get_histogram_edges`](@ref). Returns `nothing` when no positive
samples are available, otherwise returns `(edges, counts)`.
"""
function histogram_profile(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)
    edges = get_histogram_edges(samples; n_bins=n_bins)
    if edges === nothing
        return nothing
    end

    counts = zeros(Int, length(edges) - 1)
    for value in samples
        value <= 0.0 && continue
        idx = searchsortedlast(edges, value)
        idx = clamp(idx, 1, length(edges) - 1)
        if value == edges[end]
            idx = length(edges) - 1
        end
        counts[idx] += 1
    end

    return edges, counts
end

function mirrored_histogram_bars!(ax::Axis, x0::Real, samples::Vector{Float64}, side::Symbol, fill_color, edge_color;
    max_width::Float64=0.38, n_bins::Union{Nothing, Int}=nothing)
    profile = histogram_profile(samples; n_bins=n_bins)
    if profile === nothing
        return nothing
    end

    edges, counts = profile
    max_count = maximum(counts)
    if max_count == 0
        return nothing
    end

    widths = max_width .* (counts ./ max_count)
    for idx in eachindex(counts)
        widths[idx] == 0 && continue
        x1 = side == :left ? x0 - widths[idx] : x0
        x2 = side == :left ? x0 : x0 + widths[idx]
        y1 = edges[idx]
        y2 = edges[idx + 1]
        poly!(ax, Point2f[(x1, y1), (x2, y1), (x2, y2), (x1, y2)];
            color=fill_color, strokecolor=:transparent)
    end

    outer_x = Float64[]
    outer_y = Float64[]
    for idx in eachindex(counts)
        push!(outer_y, edges[idx], edges[idx + 1])
        x_outer = side == :left ? x0 - widths[idx] : x0 + widths[idx]
        push!(outer_x, x_outer, x_outer)
    end
    lines!(ax, outer_x, outer_y; color=edge_color, linewidth=2)
    return nothing
end

function add_histogram_pair!(ax::Axis, x0::Real, data_r::Vector{Float64}, data_l::Vector{Float64})
    mirrored_histogram_bars!(ax, x0, data_r, :right, PHENOM_D_FILL_COLOR, PHENOM_D_COLOR)
    mirrored_histogram_bars!(ax, x0, data_l, :left, PHENOM_HM_FILL_COLOR, PHENOM_HM_COLOR)

    combined = vcat(data_r, data_l)
    combined = combined[combined .> 0.0]
    if !isempty(combined)
        lines!(ax, [x0, x0], [minimum(combined), maximum(combined)]; color=(:black, 0.45), linewidth=1)
    end
    return nothing
end

function add_bootstrap_summary!(ax::Axis, x0::Real, samples::Vector{Float64}, side::Symbol, color)
    positive_samples = samples[samples .> 0.0]
    isempty(positive_samples) && return nothing

    q05 = quantile(positive_samples, 0.05)
    q50 = quantile(positive_samples, 0.50)
    q95 = quantile(positive_samples, 0.95)
    x_pos = side == :left ? x0 - 0.16 : x0 + 0.16

    lines!(ax, [x_pos, x_pos], [q05, q95]; color=color, linewidth=3)
    scatter!(ax, [x_pos], [q50]; color=color, markersize=16, strokecolor=:white, strokewidth=1.0)
    return nothing
end

function add_gwtc3_reference!(ax::Axis, x0::Real, pno::String)
    haskey(GWTC3_REFERENCE, pno) || return nothing
    scatter!(ax, [x0], [GWTC3_REFERENCE[pno]];
        color=GWTC3_COLOR, marker=:diamond, markersize=18,
        strokecolor=:white, strokewidth=1.0)
    return nothing
end

"""
    read_plot_data(population_results_file, config)

Read the per-waveform, per-post-Newtonian-order single-event and bootstrap
constraints from `population_results_file` using the plotting `config`.

Returns `(single_event_constraints, bootstrap_constraints)`, where each entry is
indexed as `constraints[waveform_family][pn_order]`.
"""
function read_plot_data(population_results_file::AbstractString, config::Dict)
    single_event_constraints = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_constraints = Dict{String, Dict{String, Vector{Float64}}}()

    h5open(population_results_file, "r") do file
        network_group = file[config["network"]]
        for wf_fam in config["waveform_families"]
            single_event_constraints[wf_fam] = Dict{String, Vector{Float64}}()
            bootstrap_constraints[wf_fam] = Dict{String, Vector{Float64}}()
            waveform_group = network_group[wf_fam]

            for pno in config["pn_orders"]
                pno_group = waveform_group[createSED.pnoString(pno)]
                single_event_constraints[wf_fam][pno] = read(pno_group, "single_event_constraints")
                bootstrap_constraints[wf_fam][pno] = read(pno_group, "bootstrap_constraints")
            end
        end
    end

    return single_event_constraints, bootstrap_constraints
end

function histogram_limits(
    pn_orders::AbstractVector{<:AbstractString},
    single_event_constraints::Dict{String, Dict{String, Vector{Float64}}},
    waveform_families::AbstractVector{<:AbstractString},
)
    y_min = Inf
    y_max = 0.0

    for pno in pn_orders
        for wf_fam in waveform_families
            samples = get(get(single_event_constraints, wf_fam, Dict{String, Vector{Float64}}()), pno, Float64[])
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
    single_event_constraints::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}},
    waveform_families::AbstractVector{<:AbstractString},
)
    values = Float64[]

    for pno in pn_orders
        for wf_fam in waveform_families
            single_event = get(get(single_event_constraints, wf_fam, Dict{String, Vector{Float64}}()), pno, Float64[])
            append!(values, single_event[single_event .> 0.0])

            bootstrap = get(get(bootstrap_constraints, wf_fam, Dict{String, Vector{Float64}}()), pno, Float64[])
            positive_bootstrap = bootstrap[bootstrap .> 0.0]
            if !isempty(positive_bootstrap)
                append!(values, quantile(positive_bootstrap, [0.05, 0.50, 0.95]))
            end
        end

        if haskey(GWTC3_REFERENCE, pno) && GWTC3_REFERENCE[pno] > 0.0
            push!(values, GWTC3_REFERENCE[pno])
        end
    end

    isempty(values) && return nothing
    return (minimum(values), maximum(values))
end

function expand_log_limits(y_limits::Tuple{<:Real, <:Real}; lower_pad_decades::Real=0.08, upper_pad_decades::Real=0.08)
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
    if length(exponents) > 2
        exponents = exponents[2:end-1]
    end
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
    single_event_constraints::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}})
    for (idx, pno) in enumerate(pn_orders)
        data_d = get(get(single_event_constraints, "PhenomD", Dict{String, Vector{Float64}}()), pno, Float64[])
        data_hm = get(get(single_event_constraints, "PhenomHM", Dict{String, Vector{Float64}}()), pno, Float64[])
        boot_d = get(get(bootstrap_constraints, "PhenomD", Dict{String, Vector{Float64}}()), pno, Float64[])
        boot_hm = get(get(bootstrap_constraints, "PhenomHM", Dict{String, Vector{Float64}}()), pno, Float64[])
        add_histogram_pair!(ax, idx, data_d, data_hm)
        add_bootstrap_summary!(ax, idx, boot_d, :right, PHENOM_D_COLOR)
        add_bootstrap_summary!(ax, idx, boot_hm, :left, PHENOM_HM_COLOR)
        add_gwtc3_reference!(ax, idx, pno)
    end
    return ax
end

function add_fig1_legend!(fig::Figure, target_slot)
    elements = [
        PolyElement(color=PHENOM_D_FILL_COLOR, strokecolor=PHENOM_D_COLOR),
        PolyElement(color=PHENOM_HM_FILL_COLOR, strokecolor=PHENOM_HM_COLOR),
        MarkerElement(color=:black, marker=:circle, markersize=15, strokecolor=:white, strokewidth=1.0),
        LineElement(color=:black, linewidth=3),
        MarkerElement(color=GWTC3_COLOR, marker=:diamond, markersize=18, strokecolor=:white, strokewidth=1.0),
    ]
    labels = [
        L"\text{IMRPhenomD}",
        L"\text{IMRPhenomHM}", 
        L"\text{Population constraint median}", 
        L"\text{Population constraint 90\% CI}", 
        L"\text{GWTC-3 TGR (SEOBNRv4\_ROM)}"
    ]
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

function build_figure(
    left_orders::Vector{String},
    right_orders::Vector{String},
    single_event_constraints::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}},
    aligned_left_limits::Tuple{<:Real, <:Real},
    right_scale::Real,
)
    CairoMakie.activate!()

    fig = Figure(size=FIG1_SIZE, backgroundcolor=:white)
    left_top_grid = GridLayout(fig[1, 1])
    right_top_grid = GridLayout(fig[1, 2])
    left_ax = Axis(fig[2, 1], backgroundcolor=:white, yscale=log10)
    right_ax = Axis(fig[2, 2], backgroundcolor=:white, yscale=log10)

    panel_width_weights = Float64[length(left_orders), length(right_orders)]
    width_total = sum(panel_width_weights)
    colsize!(fig.layout, 1, Relative(panel_width_weights[1] / width_total))
    colsize!(fig.layout, 2, Relative(panel_width_weights[2] / width_total))
    rowsize!(fig.layout, 1, Relative(FIG1_TOP_ROW_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG1_TOP_ROW_FRACTION))
    colgap!(fig.layout, FIG1_COL_GAP)
    rowgap!(fig.layout, FIG1_LABEL_ROW_GAP)

    build_top_label_row!(left_top_grid, left_orders)
    build_top_label_row!(right_top_grid, right_orders)

    configure_panel_axis!(left_ax, left_orders, aligned_left_limits, decade_ticks(aligned_left_limits);
        ylabel=L"|\delta\varphi_p|", show_ylabel=true)
    configure_panel_axis!(right_ax, right_orders,
        (aligned_left_limits[1] / right_scale, aligned_left_limits[2] / right_scale),
        decade_ticks((aligned_left_limits[1] / right_scale, aligned_left_limits[2] / right_scale));
        ylabel="", show_ylabel=false)

    build_panel!(left_ax, left_orders, single_event_constraints, bootstrap_constraints)
    build_panel!(right_ax, right_orders, single_event_constraints, bootstrap_constraints)
    add_fig1_legend!(fig, fig[2, 2])

    return fig
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
"""
    run_plot_fig1(config::Dict)

Generate the Figure 1 comparison plot from the population-results file
specified by `config` and save it to the configured output path.

The function reads the single-event and bootstrap constraints, builds a
dedicated panel for PN order `-1` and a second panel for the remaining PN
orders, combines both panels into the final figure, writes the figure to disk,
and returns the output file path.
"""
function run_plot_fig1(config::Dict)
    population_results_file = get_population_results_file(config)
    isfile(population_results_file) || throw(ArgumentError("Population results file does not exist: $(population_results_file)"))

    single_event_constraints, bootstrap_constraints = read_plot_data(population_results_file, config)

    left_orders = ["-1"]
    right_orders = filter(pno -> pno != "-1", config["pn_orders"])

    bootstrap_medians = Dict{String, Float64}()
    for pno in ("-1", "0")
        positive_samples = Float64[]
        for wf_fam in config["waveform_families"]
            samples = get(get(bootstrap_constraints, wf_fam, Dict{String, Vector{Float64}}()), pno, Float64[])
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

    left_hist_limits = histogram_limits(left_orders, single_event_constraints, config["waveform_families"])
    right_hist_limits = histogram_limits(right_orders, single_event_constraints, config["waveform_families"])
    left_data_limits = panel_data_limits(left_orders, single_event_constraints, bootstrap_constraints, config["waveform_families"])
    aligned_left_limits = isnothing(left_data_limits) ? (1e-5, 1.0) : (left_data_limits[1] / 10, left_data_limits[2])

    if left_hist_limits !== nothing
        aligned_left_limits = (
            min(aligned_left_limits[1], left_hist_limits[1]),
            max(aligned_left_limits[2], left_hist_limits[2]),
        )
    end

    if right_hist_limits !== nothing
        aligned_left_limits = (
            min(aligned_left_limits[1], right_hist_limits[1] * right_scale),
            max(aligned_left_limits[2], right_hist_limits[2] * right_scale),
        )
    end

    right_data_limits = panel_data_limits(right_orders, single_event_constraints, bootstrap_constraints, config["waveform_families"])
    if right_data_limits !== nothing
        aligned_left_limits = (
            min(aligned_left_limits[1], right_data_limits[1] * right_scale),
            max(aligned_left_limits[2], right_data_limits[2] * right_scale),
        )
    end

    aligned_left_limits = expand_log_limits(aligned_left_limits)

    fig = build_figure(
        left_orders,
        right_orders,
        single_event_constraints,
        bootstrap_constraints,
        aligned_left_limits,
        right_scale,
    )

    png_output_file, pdf_output_file = get_plot_output_files(config)
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
    config = read_config(config_file)

    return run_plot_fig1(config)
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
