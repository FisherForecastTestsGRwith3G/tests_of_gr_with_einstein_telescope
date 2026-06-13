using TOML
using HDF5
using CairoMakie
using LaTeXStrings
import JSON

include("_config_parser.jl")
include("_plot_style.jl")
include("../create_single_event_datasets/createSED.jl")
include("_plotting_utils.jl")

## WARNING: NOT HUMAN CONTROLLED YET! 

const TARGET_PN_ORDERS = ["-1", "1", "log(3.)"]
const TARGET_WAVEFORM_FAMILY = "PhenomHM"
const PN_COLORS = Dict(
    "-1" => FIG9_IMPROVEMENT_HIGH_COLOR,
    "1" => FIG9_IMPROVEMENT_MIDDLE_COLOR,
    "log(3.)" => FIG9_IMPROVEMENT_LOW_COLOR,
)
const WAVEFORM_LINESTYLES = Dict(
    "PhenomD" => nothing,
    "PhenomHM" => :dash,
)
const FIG_SCALING_SIZE = (1800, 900)
const FIG3_FIGURE_PADDING = (0, 40, 0, 0)
const FIG3_LEGEND_ROW_FRACTION = 0.18
const SQRT_N_GUIDE_COLOR = (:gray35, 0.65)
const CBRT_N_GUIDE_COLOR = (:gray45, 0.65)
const POPULATION_LINE_WIDTH = 3.5
const BEST_EVENT_LINE_WIDTH = 3.5
const GUIDE_LINE_WIDTH = 3
const LEGEND_POPULATION_LINE_WIDTH = 3.5
const LEGEND_BEST_EVENT_LINE_WIDTH = 3.5
const LEGEND_GUIDE_LINE_WIDTH = 3

function get_scaling_results_file(config::Dict)
    return joinpath(
        @__DIR__,
        config["scaling_outdir"],
        "scaling_results_$(config["scaling_tag"]).h5"
    )
end

function get_scaling_plot_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"], "fig_3")
    output_stem = "fig3_$(config["scaling_tag"])"
    return (
        joinpath(output_dir, output_stem * ".png"),
        joinpath(output_dir, output_stem * ".pdf"),
    )
end

function pno_plot_label(pno::AbstractString)
    pno == "-1" && return L"\varphi_{-2}"
    pno == "1" && return L"\varphi_{2}"
    pno == "log(3.)" && return L"\varphi_{6\ell}"
    return latexstring(raw"\varphi_{", pno, "}")
end

function decade_limits_and_ticks(values::Vector{Float64})
    y_min = minimum(values)
    y_max = maximum(values)
    decade_min = floor(Int, log10(y_min))
    decade_max = ceil(Int, log10(y_max))
    return (exp10(decade_min), exp10(decade_max)), decade_ticks(decade_min, decade_max)
end

"""
    read_scaling_results(scaling_results_file::AbstractString, config::Dict)

Read the Fig. 3 scaling-analysis results for `config["network"]` from the HDF5
file at `scaling_results_file`.

Returns a dictionary keyed by PN order string. Each value contains the population
sizes and the median, 5th percentile, and 95th percentile curves for the
population and best-single-detector selections. Returns `nothing` when the
target waveform family is absent.
"""
function read_scaling_results(scaling_results_file::AbstractString, config::Dict)
    results = Dict{String, Any}()

    h5open(scaling_results_file, "r") do file
        network_group = file[config["network"]]

        haskey(network_group, TARGET_WAVEFORM_FAMILY) || return nothing
        wf_group = network_group[TARGET_WAVEFORM_FAMILY]

        for pno in TARGET_PN_ORDERS
            haskey(wf_group, createSED.pnoString(pno)) || continue
            pno_group = wf_group[createSED.pnoString(pno)]

            results[pno] = (
                population_sizes=read(pno_group, "population_sizes"),
                population=(
                    median=read(pno_group["population"], "median"),
                    q05=read(pno_group["population"], "q05"),
                    q95=read(pno_group["population"], "q95"),
                ),
                best_single=(
                    median=read(pno_group["best_single"], "median"),
                    q05=read(pno_group["best_single"], "q05"),
                    q95=read(pno_group["best_single"], "q95"),
                ),
            )
        end
    end

    return results
end

"""
    plot_scaling_band!(ax::Axis, x::AbstractVector, q05::AbstractVector, q95::AbstractVector, color; alpha=0.3)

Plot the positive finite entries of the 5%-95% quantile band on `ax`.
"""
function plot_scaling_band!(ax::Axis, x::AbstractVector, q05::AbstractVector, q95::AbstractVector, color;
    alpha::Real=0.3)
    positive = positive_finite_mask(x, q05, q95)
    any(positive) || return nothing

    band!(ax, x[positive], q05[positive], q95[positive];
        color=(color, alpha),
    )
end

"""
    plot_scaling_line!(ax::Axis, x::AbstractVector, y::AbstractVector, color, linestyle; linewidth=3.0)

Plot the positive finite entries of `y` against `x` on `ax`.
"""
function plot_scaling_line!(ax::Axis, x::AbstractVector, y::AbstractVector, color, linestyle;
    linewidth::Real=3.0)
    positive = positive_finite_mask(x, y)
    any(positive) || return nothing

    lines!(ax, x[positive], y[positive];
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
    )
end

"""
    reference_power_law_curve(x::AbstractVector, y::AbstractVector, exponent::Real; extension_decades=0.18)

Return `(x_guide, y_guide, x_positive)` for a reference power-law curve anchored
at the final positive finite point in `x` and `y`, or `nothing` if no such
point exists.
"""
function reference_power_law_curve(x::AbstractVector, y::AbstractVector, exponent::Real; extension_decades::Real=0.18)
    positive = positive_finite_mask(x, y)
    any(positive) || return nothing

    x_positive = Float64.(x[positive])
    y_positive = Float64.(y[positive])
    anchor_idx = length(x_positive)
    amplitude = y_positive[anchor_idx] * x_positive[anchor_idx]^exponent
    x_min = minimum(x_positive) / exp10(extension_decades)
    x_max = maximum(x_positive) * exp10(extension_decades)
    x_guide = exp10.(range(log10(x_min), log10(x_max); length=length(x_positive) + 8))

    return x_guide, amplitude ./ (x_guide .^ exponent), x_positive
end

"""
    add_fig3_legend!(fig::Figure, target_slot)

Add the Fig. 3 PN-order and constraint-type legend to `target_slot`.

Returns `fig`.
"""
function add_fig3_legend!(fig::Figure, target_slot)
    pn_elements = [
        LineElement(color=PN_COLORS[pno], linestyle=nothing, linewidth=5)
        for pno in TARGET_PN_ORDERS
    ]
    pn_labels = Any[pno_plot_label(pno) for pno in TARGET_PN_ORDERS]
    push!(pn_elements, LineElement(color=(:transparent, 0.0), linewidth=0))
    push!(pn_labels, "")

    summary_elements = [
        [
            PolyElement(color=(:gray70, 0.28), strokecolor=:transparent),
            LineElement(color=:black, linestyle=nothing, linewidth=LEGEND_POPULATION_LINE_WIDTH),
        ],
        [
            PolyElement(color=(:gray70, 0.18), strokecolor=:transparent),
            LineElement(color=:black, linestyle=:dash, linewidth=LEGEND_BEST_EVENT_LINE_WIDTH),
        ],
        LineElement(color=SQRT_N_GUIDE_COLOR, linestyle=:dashdot, linewidth=LEGEND_GUIDE_LINE_WIDTH),
        LineElement(color=CBRT_N_GUIDE_COLOR, linestyle=:dot, linewidth=LEGEND_GUIDE_LINE_WIDTH),
    ]
    summary_labels = [
        L"\text{Population constraint}",
        L"\text{Best single event constraint}",
        L"N_{\mathrm{obs}}^{-1/2}",
        L"N_{\mathrm{obs}}^{-1/3}",
    ]

    Legend(
        target_slot,
        [pn_elements, summary_elements],
        [pn_labels, summary_labels],
        [L"\text{PN order}", L"\text{Constraint type}"];
        tellwidth=true,
        tellheight=false,
        halign=:center,
        valign=:top,
        framevisible=true,
        backgroundcolor=(:white, 0.9),
        labelsize=LEGEND_FONT_SIZE,
        titlesize=LEGEND_FONT_SIZE,
        rowgap=8,
        patchsize=(36, 24),
        colgap=18,
        margin=(22, 22, 12, 12),
        orientation=:horizontal,
    )

    return fig
end

function build_scaling_figure(results::Dict{String, Any}, config::Dict)
    fig = Figure(size=FIG_SCALING_SIZE, backgroundcolor=:white, figure_padding=FIG3_FIGURE_PADDING)
    legend_slot = fig[1, 1]
    ax = Axis(
        fig[2, 1],
        backgroundcolor=:white,
        xscale=log10,
        yscale=log10,
        xlabel=L"N_\mathrm{obs}",
        ylabel=L"|\delta\varphi_p|",
        xlabelsize=GUIDE_FONT_SIZE,
        ylabelsize=GUIDE_FONT_SIZE,
        xticklabelsize=TICK_FONT_SIZE,
        yticklabelsize=TICK_FONT_SIZE,
        spinewidth=2,
        xtickwidth=2,
        ytickwidth=2,
    )
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xminorgridvisible = true
    ax.yminorgridvisible = true
    ax.xgridcolor = (:gray70, 0.45)
    ax.ygridcolor = (:gray70, 0.45)
    ax.xminorgridcolor = (:gray70, 0.25)
    ax.yminorgridcolor = (:gray70, 0.25)
    ax.xminorticks = IntervalsBetween(9)
    ax.yminorticks = IntervalsBetween(9)
    ax.xtickalign = 1
    ax.ytickalign = 1
    colsize!(fig.layout, 1, Relative(1.0))
    rowsize!(fig.layout, 1, Relative(FIG3_LEGEND_ROW_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG3_LEGEND_ROW_FRACTION))
    rowgap!(fig.layout, 0.8)

    plotted_any = false
    sqrt_guide_label_added = false
    cbrt_guide_label_added = false
    x_min = Inf
    x_max = -Inf
    y_values = Float64[]

    for pno in TARGET_PN_ORDERS
        haskey(results, pno) || continue
        result = results[pno]
        color = PN_COLORS[pno]

        plot_scaling_band!(
            ax,
            result.population_sizes,
            result.population.q05,
            result.population.q95,
            color,
        )
        plot_scaling_line!(
            ax,
            result.population_sizes,
            result.population.median,
            color,
            nothing,
            linewidth=POPULATION_LINE_WIDTH,
        )
        plot_scaling_band!(
            ax,
            result.population_sizes,
            result.best_single.q05,
            result.best_single.q95,
            color;
            alpha=0.15,
        )
        plot_scaling_line!(
            ax,
            result.population_sizes,
            result.best_single.median,
            color,
            :dash,
            linewidth=BEST_EVENT_LINE_WIDTH,
        )

        sqrt_guide_curve = reference_power_law_curve(result.population_sizes, result.population.median, 0.5)
        if sqrt_guide_curve !== nothing
            guide_x, guide_y, guide_data_x = sqrt_guide_curve
            x_min = min(x_min, minimum(guide_data_x))
            x_max = max(x_max, maximum(guide_data_x))
            lines!(ax, guide_x, guide_y;
                color=SQRT_N_GUIDE_COLOR,
                linestyle=:dashdot,
                linewidth=GUIDE_LINE_WIDTH,
            )
            if !sqrt_guide_label_added
                text!(
                    ax,
                    guide_x[end],
                    guide_y[end];
                    text=L"N_{\mathrm{obs}}^{-1/2}",
                    fontsize=LEGEND_FONT_SIZE,
                    color=SQRT_N_GUIDE_COLOR,
                    align=(:right, :bottom),
                    offset=(-8, 8),
                )
                sqrt_guide_label_added = true
            end
        end

        cbrt_guide_curve = reference_power_law_curve(result.population_sizes, result.best_single.median, 1 / 3)
        if cbrt_guide_curve !== nothing
            guide_x, guide_y, guide_data_x = cbrt_guide_curve
            x_min = min(x_min, minimum(guide_data_x))
            x_max = max(x_max, maximum(guide_data_x))
            lines!(ax, guide_x, guide_y;
                color=CBRT_N_GUIDE_COLOR,
                linestyle=:dot,
                linewidth=GUIDE_LINE_WIDTH,
            )
            if !cbrt_guide_label_added
                text!(
                    ax,
                    guide_x[end],
                    guide_y[end];
                    text=L"N_{\mathrm{obs}}^{-1/3}",
                    fontsize=LEGEND_FONT_SIZE,
                    color=CBRT_N_GUIDE_COLOR,
                    align=(:right, :bottom),
                    offset=(-8, 8),
                )
                cbrt_guide_label_added = true
            end
        end

        append_positive_finite_values!(y_values, result.population_sizes, result.population.q05)
        append_positive_finite_values!(y_values, result.population_sizes, result.population.q95)
        append_positive_finite_values!(y_values, result.population_sizes, result.best_single.q05)
        append_positive_finite_values!(y_values, result.population_sizes, result.best_single.median)
        append_positive_finite_values!(y_values, result.population_sizes, result.best_single.q95)
        plotted_any = true
    end

    plotted_any || throw(ArgumentError("No scaling results found for target PN orders: $(join(TARGET_PN_ORDERS, ", "))"))

    add_fig3_legend!(fig, legend_slot)

    if isfinite(x_min) && isfinite(x_max)
        xlims!(ax, x_min, x_max)
        ax.xticks = decade_ticks((x_min, x_max); trim_outer=false)
    end

    if !isempty(y_values)
        y_limits, y_tick_spec = decade_limits_and_ticks(y_values)
        ylims!(ax, y_limits...)
        ax.yticks = y_tick_spec
    end

    return fig
end

function run_plot_scaling_analysis(config::Dict)
    scaling_results_file = get_scaling_results_file(config)
    isfile(scaling_results_file) || throw(ArgumentError("Scaling results file does not exist: $(scaling_results_file)"))

    results = read_scaling_results(scaling_results_file, config)
    fig = build_scaling_figure(results, config)

    png_output_file, pdf_output_file = get_scaling_plot_output_files(config)
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

    return run_plot_scaling_analysis(config)
end

fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))


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
