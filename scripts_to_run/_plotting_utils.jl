const GWTC3_COLOR = :black
const GUIDE_FONT_SIZE = 32
const TICK_FONT_SIZE = 32
const LEGEND_FONT_SIZE = 28
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

"""
    get_histogram_edges(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)

Return logarithmically spaced histogram bin edges for the positive entries in
`samples`, or `nothing` when no positive samples are available.
"""
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

"""
    histogram_profile(samples::Vector{Float64}; n_bins::Union{Nothing, Int}=nothing)

Build histogram bin counts using logarithmically spaced positive-sample edges.
Returns `nothing` when no positive samples are available.
"""
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

function add_gwtc3_reference!(ax::Axis, x0::Real, pno::String)
    haskey(GWTC3_REFERENCE, pno) || return nothing
    scatter!(ax, [x0], [GWTC3_REFERENCE[pno]];
        color=GWTC3_COLOR, marker=:diamond, markersize=18,
        strokecolor=:white, strokewidth=1.0)
    return nothing
end

function histogram_limits(
    pn_orders::AbstractVector{<:AbstractString},
    constraints_by_group::Dict{String, Dict{String, Vector{Float64}}},
    group_names::AbstractVector{<:AbstractString},
)
    y_min = Inf
    y_max = 0.0

    for pno in pn_orders
        for group_name in group_names
            samples = get(get(constraints_by_group, group_name, Dict{String, Vector{Float64}}()), pno, Float64[])
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
    single_event_by_group::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_by_group::Dict{String, Dict{String, Vector{Float64}}},
    group_names::AbstractVector{<:AbstractString},
)
    values = Float64[]

    for pno in pn_orders
        for group_name in group_names
            single_event = get(get(single_event_by_group, group_name, Dict{String, Vector{Float64}}()), pno, Float64[])
            append!(values, single_event[single_event .> 0.0])

            bootstrap = get(get(bootstrap_by_group, group_name, Dict{String, Vector{Float64}}()), pno, Float64[])
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
