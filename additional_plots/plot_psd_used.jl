using CairoMakie
using ColorSchemes
using DelimitedFiles
using LaTeXStrings

const PSD_DATA_DIR = joinpath(@__DIR__, "..", "create_single_event_datasets", "psd_data")
const HLV_DATA_DIR = joinpath(PSD_DATA_DIR, "hlv_curves", "O3b")
const ET_DATA_DIR = joinpath(PSD_DATA_DIR, "et_curves")
const OUTPUT_DIR = joinpath(@__DIR__, "results")

const FIGURE_SIZE = (1000, 500)
const X_LIMITS = (0.8, 12000)
const Y_LIMITS = (1e-25, 1e-16)
const ET_FMIN = 2
const LVK_FMIN = 20
const LVK_FMAX = 5000
const AXIS_LABEL_SIZE = 24
const TICK_LABEL_SIZE = 18
const TICK_SIZE = 6
const CURVE_LINE_WIDTH = 1.5
const LEGEND_LINE_WIDTH = 3 * CURVE_LINE_WIDTH
const CUTOFF_LINE_COLOR = (:gray45, 0.8)
const CUTOFF_LINE_STYLE = :dash
const CUTOFF_LINE_WIDTH = 1.5

const COMMON_PALETTE = ColorSchemes.seaborn_colorblind
const CURVE_COLORS = Dict(
    :L1 => COMMON_PALETTE[1],
    :H1 => COMMON_PALETTE[2],
    :V1 => COMMON_PALETTE[3],
    :ET10 => COMMON_PALETTE[4],
    :ET15 => COMMON_PALETTE[5],
)

const CURVE_LABELS = Dict(
    :L1 => L"\mathrm{LIGO\,O3b\,L}",
    :H1 => L"\mathrm{LIGO\,O3b\,H}",
    :V1 => L"\mathrm{Virgo\,O3b\,V}",
    :ET10 => L"\mathrm{ET}\ 10\mathrm{km}",
    :ET15 => L"\mathrm{ET}\ 15\mathrm{km}",
)

const ASD_FILES = Dict(
    :V1 => joinpath(HLV_DATA_DIR, "O3-V1-1265246178_sensitivity_strain_asd.txt"),
    :H1 => joinpath(HLV_DATA_DIR, "O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt"),
    :L1 => joinpath(HLV_DATA_DIR, "O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt"),
)

const PSD_FILES = Dict(
    :ET10 => joinpath(ET_DATA_DIR, "ET10km.txt"),
    :ET15 => joinpath(ET_DATA_DIR, "ET15km.txt"),
)

function loadtxt(path::AbstractString)
    return Matrix{Float64}(readdlm(path, Float64))
end

function plot_curve_segments!(ax::Axis, data::AbstractMatrix{<:Real}, selected::AbstractVector{Bool};
    label=nothing, color, linewidth::Real=CURVE_LINE_WIDTH, faded_alpha::Real=0.4, transform_y=identity)

    x = view(data, :, 1)
    y = transform_y.(view(data, :, 2))

    lines!(
        ax,
        x[selected],
        y[selected];
        color=(color, 1.0),
        linewidth=linewidth,
        label=label,
    )
    lines!(
        ax,
        x[.!selected],
        y[.!selected];
        color=(color, faded_alpha),
        linewidth=linewidth,
    )

    return nothing
end

function add_curve_legend!(target_slot)
    curve_order = (:L1, :H1, :V1, :ET10, :ET15)
    elements = [
        LineElement(color=CURVE_COLORS[curve], linestyle=nothing, linewidth=LEGEND_LINE_WIDTH)
        for curve in curve_order
    ]
    labels = [CURVE_LABELS[curve] for curve in curve_order]

    return Legend(
        target_slot,
        elements,
        labels;
        tellwidth=false,
        tellheight=false,
        halign=:right,
        valign=:top,
        framevisible=true,
        labelsize=22,
        patchsize=(48, 24),
        margin=(10, 10, 10, 10),
    )
end

function build_plot()
    CairoMakie.activate!()

    asd_data = Dict(detector => loadtxt(path) for (detector, path) in ASD_FILES)
    psd_data = Dict(detector => loadtxt(path) for (detector, path) in PSD_FILES)

    fig = Figure(size=FIGURE_SIZE, backgroundcolor=:white)
    ax = Axis(
        fig[1, 1];
        xscale=log10,
        yscale=log10,
        xlabel=L"\mathrm{Frequency}\,\mathrm{[Hz]}",
        ylabel=L"\mathrm{Strain}\,[1/\sqrt{\mathrm{Hz}}]",
        xlabelsize=AXIS_LABEL_SIZE,
        ylabelsize=AXIS_LABEL_SIZE,
        xticklabelsize=TICK_LABEL_SIZE,
        yticklabelsize=TICK_LABEL_SIZE,
        xticksize=TICK_SIZE,
        yticksize=TICK_SIZE,
        xgridvisible=true,
        ygridvisible=true,
        xminorgridvisible=true,
        yminorgridvisible=true,
        xgridwidth=0.2,
        ygridwidth=0.2,
        xminorgridwidth=0.2,
        yminorgridwidth=0.2,
        xminorticks=IntervalsBetween(9),
        yminorticks=IntervalsBetween(9),
        xminorticksvisible=true,
        yminorticksvisible=true,
        xticksvisible=true,
        yticksvisible=true,
        xtickalign=1,
        ytickalign=1,
    )

    plot_curve_segments!(
        ax,
        asd_data[:L1],
        asd_data[:L1][:, 1] .> LVK_FMIN;
        label=CURVE_LABELS[:L1],
        color=CURVE_COLORS[:L1],
    )
    plot_curve_segments!(
        ax,
        asd_data[:H1],
        asd_data[:H1][:, 1] .> LVK_FMIN;
        label=CURVE_LABELS[:H1],
        color=CURVE_COLORS[:H1],
    )
    plot_curve_segments!(
        ax,
        asd_data[:V1],
        (asd_data[:V1][:, 1] .> LVK_FMIN) .& (asd_data[:V1][:, 1] .< LVK_FMAX);
        label=CURVE_LABELS[:V1],
        color=CURVE_COLORS[:V1],
    )
    plot_curve_segments!(
        ax,
        psd_data[:ET10],
        psd_data[:ET10][:, 1] .> ET_FMIN;
        label=CURVE_LABELS[:ET10],
        color=CURVE_COLORS[:ET10],
        transform_y=sqrt,
    )
    plot_curve_segments!(
        ax,
        psd_data[:ET15],
        psd_data[:ET15][:, 1] .> ET_FMIN;
        label=CURVE_LABELS[:ET15],
        color=CURVE_COLORS[:ET15],
        transform_y=sqrt,
    )

    vlines!(
        ax,
        [ET_FMIN, LVK_FMIN, LVK_FMAX];
        color=CUTOFF_LINE_COLOR,
        linestyle=CUTOFF_LINE_STYLE,
        linewidth=CUTOFF_LINE_WIDTH,
    )

    xlims!(ax, X_LIMITS...)
    ylims!(ax, Y_LIMITS...)
    add_curve_legend!(fig[1, 1])

    return fig
end

function main()
    mkpath(OUTPUT_DIR)
    fig = build_plot()
    save(joinpath(OUTPUT_DIR, "psds_used.pdf"), fig)
    return nothing
end

main()
