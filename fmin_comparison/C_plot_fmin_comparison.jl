using HDF5
using CairoMakie
using LaTeXStrings
using Statistics

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

const ERROR_BARS_INTERVAL_MULTIPLIER = 1.645
const FMIN_COLORS = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7"]
const FMIN_MARKERS = [:circle, :rect, :diamond, :utriangle, :dtriangle]
const FIGURE_SIZE = (1800, 760)
const FONT_SIZE_LABEL = 32
const FONT_SIZE_TICK = 28
const FONT_SIZE_LEGEND = 25

const PN_LABELS = Dict(
    "-1" => L"\varphi_{-2}",
    "0" => L"\varphi_{0}",
    "0.5" => L"\varphi_{1}",
    "1" => L"\varphi_{2}",
    "1.5" => L"\varphi_{3}",
    "2" => L"\varphi_{4}",
    "log(2.5)" => L"\varphi_{5\,\ell}",
    "3" => L"\varphi_{6}",
    "log(3.)" => L"\varphi_{6\,\ell}",
    "3.5" => L"\varphi_{7}",
)

const ET_BLUEBOOK_REFERENCE = Dict(
    "-1" => 4.5e-5,
    "0" => 0.0035,
    "0.5" => 0.011,
    "1" => 0.008,
    "1.5" => 0.004,
    "2" => 0.034,
    "log(2.5)" => 0.012,
    "3" => 0.016,
    "log(3.)" => 0.088,
    "3.5" => 0.048,
)

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_plot_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"], "fmin_comparison")
    stem = "fmin_comparison_$(config["plot_tag"])"
    return (
        joinpath(output_dir, stem * ".png"),
        joinpath(output_dir, stem * ".pdf"),
    )
end

function pretty_fmin(fmin::Real)
    value = Float64(fmin)
    if isapprox(value, round(value); atol=1e-9)
        return string(round(Int, value))
    end
    return string(value)
end

function fmin_legend_label(fmin::Real)
    return latexstring("f_{\\min} = $(pretty_fmin(fmin))\\,\\mathrm{Hz}")
end

function pn_label(pno::AbstractString)
    return get(PN_LABELS, pno, latexstring(pno))
end

function relative_x_offset(index::Integer, total::Integer)
    total == 1 && return 0.0
    return ((index - 1) / (total - 1)) * 2.0 - 1.0
end

function waveform_group(file, config::Dict, pno::AbstractString)
    pno_group_name = createSED.pnoString(pno)
    return file[pno_group_name][config["network"]][config["plot_waveform"]]
end

function shared_selection(file, config::Dict)
    selection = nothing

    for pno in config["pn_orders"]
        wf_grp = waveform_group(file, config, pno)
        snr = Float64.(read(wf_grp, "snr"))
        isnr = Float64.(read(wf_grp, "isnr"))
        invc = Bool.(read(wf_grp, "invc"))
        delta_k = Float64.(read(wf_grp, "delta_k"))

        pno_selection = BitVector(
            (snr .> config["snr_threshold"]) .&
            (isnr .> config["snr_inspiral_threshold"]) .&
            invc .&
            isfinite.(delta_k) .&
            (delta_k .> 0.0)
        )

        selection = isnothing(selection) ? pno_selection : BitVector(selection .& pno_selection)
    end

    return selection
end

function conditioned_limits_by_realization(
    dphi_k::Vector{Float64},
    delta_k::Vector{Float64},
    selection::BitVector,
    events_per_realization::Integer,
)
    selected_dphi = dphi_k[selection]
    selected_delta = delta_k[selection]
    valid = isfinite.(selected_dphi) .& isfinite.(selected_delta) .& (selected_delta .> 0.0)
    selected_dphi = selected_dphi[valid]
    selected_delta = selected_delta[valid]

    n_realizations = fld(length(selected_dphi), events_per_realization)
    n_realizations > 0 ||
        throw(ArgumentError("Not enough selected events to build one conditioned realization."))

    limits = Vector{Float64}(undef, n_realizations)
    for idx in 1:n_realizations
        first_idx = (idx - 1) * events_per_realization + 1
        last_idx = idx * events_per_realization
        dphi_chunk = selected_dphi[first_idx:last_idx]
        delta_chunk = selected_delta[first_idx:last_idx]
        mu_eff, sigma_eff = HierDist.muStdEff4ProdNormal(dphi_chunk, delta_chunk)
        limits[idx] = HierDist.estimate0Symmetric90CiGaussian(mu_eff, sigma_eff)
    end

    return limits
end

function summarize_limits(limits::Vector{Float64}; log_space::Bool=true)
    clean_limits = limits[isfinite.(limits) .& (limits .> 0.0)]
    isempty(clean_limits) && throw(ArgumentError("No positive finite upper limits available."))

    if length(clean_limits) == 1
        return (center=only(clean_limits), lower=0.0, upper=0.0, n_realizations=1)
    end

    if log_space
        log_limits = log.(clean_limits)
        center = exp(mean(log_limits))
        log_sigma = std(log_limits)
        lower = center - exp(log(center) - ERROR_BARS_INTERVAL_MULTIPLIER * log_sigma)
        upper = exp(log(center) + ERROR_BARS_INTERVAL_MULTIPLIER * log_sigma) - center
        return (center=center, lower=lower, upper=upper, n_realizations=length(clean_limits))
    end

    center = mean(clean_limits)
    sigma = std(clean_limits)
    return (
        center=center,
        lower=ERROR_BARS_INTERVAL_MULTIPLIER * sigma,
        upper=ERROR_BARS_INTERVAL_MULTIPLIER * sigma,
        n_realizations=length(clean_limits),
    )
end

function read_fmin_summaries(base_config::Dict)
    summaries = Dict{String, Dict{String, NamedTuple}}()
    selected_counts = Dict{String, Int}()

    for config in iter_fmin_configs(base_config)
        fisher_results_file = get_fisher_results_file(config)
        isfile(fisher_results_file) ||
            throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

        label = config["fmin_label"]
        summaries[label] = Dict{String, NamedTuple}()

        h5open(fisher_results_file, "r") do file
            selection = shared_selection(file, config)
            selected_counts[label] = sum(selection)

            for pno in config["pn_orders"]
                wf_grp = waveform_group(file, config, pno)
                dphi_k = Float64.(read(wf_grp, "dphi_k"))
                delta_k = Float64.(read(wf_grp, "delta_k"))
                limits = conditioned_limits_by_realization(
                    dphi_k,
                    delta_k,
                    selection,
                    config["plot_events_per_realization"],
                )
                summaries[label][pno] = summarize_limits(
                    limits;
                    log_space=config["plot_compute_log_summary"],
                )
            end
        end
    end

    return summaries, selected_counts
end

function panel_values(summaries, fmin_configs, pn_orders)
    values = Float64[]

    for config in fmin_configs
        label = config["fmin_label"]
        for pno in pn_orders
            summary = summaries[label][pno]
            push!(values, summary.center)
            summary.lower > 0.0 && push!(values, max(summary.center - summary.lower, eps()))
            summary.upper > 0.0 && push!(values, summary.center + summary.upper)
        end
    end

    for pno in pn_orders
        if haskey(ET_BLUEBOOK_REFERENCE, pno)
            push!(values, ET_BLUEBOOK_REFERENCE[pno])
        end
    end

    return values[isfinite.(values) .& (values .> 0.0)]
end

function expand_log_limits(values::Vector{Float64}; lower_pad_decades=0.08, upper_pad_decades=0.08)
    isempty(values) && return (1e-6, 1.0)
    log_min = log10(minimum(values))
    log_max = log10(maximum(values))
    return (exp10(log_min - lower_pad_decades), exp10(log_max + upper_pad_decades))
end

function configure_axis!(ax::Axis, pn_orders::Vector{String}, y_limits; ylabel::Union{String, LaTeXString}="")
    ax.xlabel = "PN order"
    ax.ylabel = ylabel
    ax.xlabelsize = FONT_SIZE_LABEL
    ax.ylabelsize = FONT_SIZE_LABEL
    ax.xticks = (collect(1:length(pn_orders)), pn_label.(pn_orders))
    ax.xticklabelsize = FONT_SIZE_TICK
    ax.yticklabelsize = FONT_SIZE_TICK
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.yminorgridvisible = true
    ax.yminorgridcolor = (:gray70, 0.35)
    ax.yminorticks = IntervalsBetween(9)
    ax.xtickalign = 1
    ax.ytickalign = 1
    xlims!(ax, 0.5, length(pn_orders) + 0.5)
    ylims!(ax, y_limits...)
    return ax
end

function plot_panel!(ax::Axis, summaries, fmin_configs, pn_orders; add_labels::Bool=false)
    offset_scale = 0.25
    n_fmin = length(fmin_configs)

    for (idx_config, config) in enumerate(fmin_configs)
        label = config["fmin_label"]
        color = FMIN_COLORS[mod1(idx_config, length(FMIN_COLORS))]
        marker = FMIN_MARKERS[mod1(idx_config, length(FMIN_MARKERS))]
        offset = offset_scale * relative_x_offset(idx_config, n_fmin)
        x_values = collect(1:length(pn_orders)) .+ offset
        y_values = [summaries[label][pno].center for pno in pn_orders]
        lower = [summaries[label][pno].lower for pno in pn_orders]
        upper = [summaries[label][pno].upper for pno in pn_orders]

        errorbars!(ax, x_values, y_values, lower, upper;
            color=color, whiskerwidth=12, linewidth=2)
        if add_labels
            scatter!(ax, x_values, y_values;
                color=color,
                marker=marker,
                markersize=18,
                strokecolor=:black,
                strokewidth=1,
                label=fmin_legend_label(config["fmin"]))
        else
            scatter!(ax, x_values, y_values;
                color=color,
                marker=marker,
                markersize=18,
                strokecolor=:black,
                strokewidth=1)
        end
    end

    reference_values = [ET_BLUEBOOK_REFERENCE[pno] for pno in pn_orders if haskey(ET_BLUEBOOK_REFERENCE, pno)]
    reference_x = [idx for (idx, pno) in enumerate(pn_orders) if haskey(ET_BLUEBOOK_REFERENCE, pno)]
    if !isempty(reference_values)
        if add_labels
            scatter!(ax, reference_x, reference_values;
                color=:black,
                marker=:diamond,
                markersize=20,
                strokecolor=:white,
                strokewidth=1,
                label=L"GW150914-like, ET $\Delta$ injection")
        else
            scatter!(ax, reference_x, reference_values;
                color=:black,
                marker=:diamond,
                markersize=20,
                strokecolor=:white,
                strokewidth=1)
        end
    end

    return ax
end

function build_figure(base_config::Dict, summaries, selected_counts)
    CairoMakie.activate!()

    fmin_configs = iter_fmin_configs(base_config)
    left_orders = ["-1"]
    right_orders = filter(pno -> pno != "-1", base_config["pn_orders"])

    left_limits = isnothing(base_config["plot_y_limits_minus_one"]) ?
        expand_log_limits(panel_values(summaries, fmin_configs, left_orders)) :
        base_config["plot_y_limits_minus_one"]
    right_limits = isnothing(base_config["plot_y_limits_remaining"]) ?
        expand_log_limits(panel_values(summaries, fmin_configs, right_orders)) :
        base_config["plot_y_limits_remaining"]

    fig = Figure(size=FIGURE_SIZE, backgroundcolor=:white)
    left_ax = Axis(fig[1, 1], backgroundcolor=:white, yscale=log10)
    right_ax = Axis(fig[1, 2], backgroundcolor=:white, yscale=log10)

    colsize!(fig.layout, 1, Relative(0.125))
    colsize!(fig.layout, 2, Relative(0.875))
    colgap!(fig.layout, 20)

    configure_axis!(left_ax, left_orders, left_limits; ylabel=L"|\delta\varphi_{\!p}|")
    configure_axis!(right_ax, right_orders, right_limits; ylabel="")

    plot_panel!(left_ax, summaries, fmin_configs, left_orders; add_labels=false)
    plot_panel!(right_ax, summaries, fmin_configs, right_orders; add_labels=true)

    axislegend(right_ax;
        position=:rb,
        framevisible=true,
        backgroundcolor=(:white, 0.92),
        labelsize=FONT_SIZE_LEGEND)

    for config in fmin_configs
        println("Selected $(selected_counts[config["fmin_label"]]) events for fmin=$(config["fmin"]) Hz")
    end

    return fig
end

function run_plot_fmin_comparison(base_config::Dict)
    summaries, selected_counts = read_fmin_summaries(base_config)
    fig = build_figure(base_config, summaries, selected_counts)

    png_output_file, pdf_output_file = get_plot_output_files(base_config)
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
    base_config = read_fmin_config(config_file)

    return run_plot_fmin_comparison(base_config)
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
