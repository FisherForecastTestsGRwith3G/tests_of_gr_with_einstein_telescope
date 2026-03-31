using TOML
using HDF5
using Plots
using LaTeXStrings
using Statistics
import JSON

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const PHENOM_D_COLOR  = "#ff7f0e"
const PHENOM_HM_COLOR = "#1f77b4"
const GWTC3_COLOR     = :black
const TITLE_FONT_SIZE = 28
const GUIDE_FONT_SIZE = 24
const TICK_FONT_SIZE  = 22
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
        return exp10.(range(log10(y_min) - 0.25, log10(y_max) + 0.25, length = 11))
    end

    if isnothing(n_bins)
        n_bins = clamp(round(Int, sqrt(length(positive_samples))), 8, 40)
    end

    return exp10.(range(log10(y_min), log10(y_max), length = n_bins + 1))
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

"""
    mirrored_histogram_bars!(plt, x0, samples, side, color; max_width=0.38, n_bins=nothing)

Draw a mirrored horizontal histogram of `samples` onto `plt`.

The histogram is anchored at the vertical line `x = x0` and drawn either to the
left or right depending on `side`, which should be `:left` or `:right`. Each bar
spans one histogram bin in the y-direction and is filled with `color`. The bar
widths are normalized by the largest bin count so that the widest bar has width
`max_width`.

Arguments:
- `plt`: plot object that is updated in place with the histogram bars.
- `x0`: x-position of the center line from which the mirrored bars extend.
- `samples`: one-dimensional sample values used to build the histogram along the
  y-axis. The function assumes these values are finite and non-negative; they are
  binned into contiguous intervals, and the resulting counts determine the bar
  widths after normalization. An empty vector is allowed and leaves `plt`
  unchanged.
- `side`: selects the drawing direction, expected to be `:left` or `:right`.
- `color`: fill color used for every histogram bar.

Keyword arguments:
- `max_width`: maximum horizontal extent of the widest histogram bar.
- `n_bins`: optional number of histogram bins; if omitted, an automatic choice is used.

Returns `nothing` when `samples` is empty or all bins have zero count; otherwise
the bars are added to `plt` in place.
"""
function mirrored_histogram_bars!(plt, x0::Real, samples::Vector{Float64}, side::Symbol, color; max_width::Float64 = 0.38, n_bins::Union{Nothing, Int}=nothing)
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
        rect = Shape([x1, x2, x2, x1], [y1, y1, y2, y2])
        plot!(plt, rect, c=color, linecolor=false, fillalpha=0.45, label=false)
    end

    outer_x = Float64[]
    outer_y = Float64[]
    for idx in eachindex(counts)
        push!(outer_y, edges[idx], edges[idx + 1])
        x_outer = side == :left ? x0 - widths[idx] : x0 + widths[idx]
        push!(outer_x, x_outer, x_outer)
    end
    plot!(plt, outer_x, outer_y, color=color, lw=2.0, label=false)
end

"""
    add_histogram_pair!(plt, x0, data_d, data_l)

Add mirrored histograms for the `data_r` and `data_l` samples at `x0`, with a
faint vertical line spanning their combined positive range.
"""
function add_histogram_pair!(plt, x0::Real, data_r::Vector{Float64}, data_l::Vector{Float64})
    mirrored_histogram_bars!(plt, x0, data_r, :right, PHENOM_D_COLOR)
    mirrored_histogram_bars!(plt, x0, data_l, :left, PHENOM_HM_COLOR)

    combined = vcat(data_r, data_l)
    combined = combined[combined .> 0.0]
    if !isempty(combined)
        plot!(plt, [x0, x0], [minimum(combined), maximum(combined)], color=:black, lw=1.0, alpha=0.45, label=false)
    end
end

"""
    add_bootstrap_summary!(plt, x0, samples, side, color)

Add a compact bootstrap summary for `samples` at horizontal position `x0` on
`plt`, drawing the 5th to 95th percentile interval and the median on the
requested `side`.
"""
function add_bootstrap_summary!(plt, x0::Real, samples::Vector{Float64}, side::Symbol, color)
    positive_samples = samples[samples .> 0.0]
    isempty(positive_samples) && return

    q05 = quantile(positive_samples, 0.05)
    q50 = quantile(positive_samples, 0.50)
    q95 = quantile(positive_samples, 0.95)
    x_pos = side == :left ? x0 - 0.16 : x0 + 0.16

    plot!(plt, [x_pos, x_pos], [q05, q95], color=color, lw=3.0, alpha=0.95, label=false)
    scatter!(plt, [x_pos], [q50], color=color, markerstrokecolor=:white, markerstrokewidth=1.0, markersize=7.875, label=false)
end

"""
    add_gwtc3_reference!(plt, x0, pno)

Add the GWTC-3 reference value for post-Newtonian order `pno` at horizontal
position `x0` on `plt`, if a reference entry is available.
"""
function add_gwtc3_reference!(plt, x0::Real, pno::String)
    haskey(GWTC3_REFERENCE, pno) || return
    scatter!(plt, [x0], [GWTC3_REFERENCE[pno]], color=GWTC3_COLOR, markercolor=GWTC3_COLOR, markerstrokecolor=:white, markerstrokewidth=0.8, markershape=:diamond, markersize=9.28125, label=false)
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

"""
    make_panel(pn_orders, single_event_constraints, bootstrap_constraints; title_text, width_scale)

Build a plot panel for the requested post-Newtonian orders using the
single-event and bootstrap constraints for each waveform family.

For each entry in `pn_orders`, the function adds mirrored histograms for the
`"PhenomD"` and `"PhenomHM"` single-event samples, overlays the corresponding
bootstrap summaries, and places the GWTC-3 reference marker when available.
The constraint dictionaries are expected to be indexed as
`constraints[waveform_family][pn_order]`.

Returns the configured plot object.
"""
function make_panel(
    pn_orders::Vector{String},
    single_event_constraints::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}};
    title_text::String,
    width_scale::Real,
)
    xticks = collect(1:length(pn_orders))
    labels = createSED.pnoLatex.(pn_orders)

    plt = plot(
        legend=false,
        yscale=:log10,
        xlabel="PN order",
        ylabel="",
        title="",
        grid=true,
        minorgrid=true,
        framestyle=:box,
        xlim=(0.5, length(pn_orders) + 0.5),
        xticks=(xticks, labels),
        size=(round(Int, 900 * width_scale), 760),
        dpi=200,
        titlefont=font(TITLE_FONT_SIZE),
        guidefont=font(GUIDE_FONT_SIZE),
        tickfont=font(TICK_FONT_SIZE),
        bottom_margin=18Plots.mm,
        top_margin=8Plots.mm,
        left_margin=8Plots.mm,
        right_margin=8Plots.mm,
    )

    for (idx, pno) in enumerate(pn_orders)
        data_d = get(get(single_event_constraints, "PhenomD", Dict{String, Vector{Float64}}()), pno, Float64[])
        data_hm = get(get(single_event_constraints, "PhenomHM", Dict{String, Vector{Float64}}()), pno, Float64[])
        boot_d = get(get(bootstrap_constraints, "PhenomD", Dict{String, Vector{Float64}}()), pno, Float64[])
        boot_hm = get(get(bootstrap_constraints, "PhenomHM", Dict{String, Vector{Float64}}()), pno, Float64[])
        add_histogram_pair!(plt, idx, data_d, data_hm)
        add_bootstrap_summary!(plt, idx, boot_d, :right, PHENOM_D_COLOR)
        add_bootstrap_summary!(plt, idx, boot_hm, :left, PHENOM_HM_COLOR)
        add_gwtc3_reference!(plt, idx, pno)
    end

    return plt
end

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

    left_panel = make_panel(left_orders, single_event_constraints, bootstrap_constraints; title_text="PN order -1", width_scale=0.45)
    right_panel = make_panel(right_orders, single_event_constraints, bootstrap_constraints; title_text="Higher PN orders", width_scale=1.35)
    left_limits = ylims(left_panel)
    aligned_left_limits = (left_limits[1] / 10, left_limits[2])

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

    ylims!(left_panel, aligned_left_limits...)
    ylims!(right_panel, aligned_left_limits[1] / right_scale, aligned_left_limits[2] / right_scale)

    final_plot = plot(
        left_panel,
        right_panel;
        layout=grid(1, 2, widths=[0.18, 0.82]),
        size=(1800, 760),
        bottom_margin=18Plots.mm,
        top_margin=8Plots.mm,
        left_margin=8Plots.mm,
        right_margin=8Plots.mm,
    )

    plot!(final_plot[2], [NaN], [NaN], seriestype=:shape, c=PHENOM_D_COLOR, linecolor=false, fillalpha=0.45, label="PhenomD")
    plot!(final_plot[2], [NaN], [NaN], seriestype=:shape, c=PHENOM_HM_COLOR, linecolor=false, fillalpha=0.45, label="PhenomHM")
    scatter!(final_plot[2], [NaN], [NaN], color=:black, markersize=7.875, label="Bootstrap median")
    plot!(final_plot[2], [NaN, NaN], [NaN, NaN], color=:black, lw=3.0, label="Bootstrap 90% CI")
    scatter!(final_plot[2], [NaN], [NaN], color=GWTC3_COLOR, markercolor=GWTC3_COLOR, markerstrokecolor=:white, markerstrokewidth=0.8, markershape=:diamond, markersize=9.28125, label="GWTC-3 TGR")
    plot!(final_plot[2], legend=:bottomright, legendfontsize=18)

    mkpath(dirname(get_plot_output_file(config)))
    output_file = get_plot_output_file(config)
    savefig(final_plot, output_file)
    println("Saved figure to $(output_file)")

    return output_file
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

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
