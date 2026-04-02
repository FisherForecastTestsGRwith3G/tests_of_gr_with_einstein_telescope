using TOML
using HDF5
using KernelDensity
using Plots
using LaTeXStrings
using Colors

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const CONTOUR_FILL_COLORS = [:aliceblue, :lightblue, :cornflowerblue]
const CONTOUR_LINE_COLOR = :royalblue4
const OBSERVABLE_COLORMAP = :viridis
const IMPROVEMENT_COLORBAR_TITLE = L"\log_{10}\!\left(\Delta k_{\mathrm{HM}} / \Delta k_{\mathrm{D}}\right)"
const FIG9_Z_THRESHOLD = 0.5
const FIG9_MC_LIMITS = (5.0, 80.0)

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_plot_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"], "fig_9")
    stem = "fig9_$(config["plot_tag"])"
    return (
        joinpath(output_dir, stem * ".png"),
        joinpath(output_dir, stem * ".pdf"),
    )
end

function get_target_pno(config::Dict)
    if "3.5" in config["pn_orders"]
        return "3.5"
    elseif length(config["pn_orders"]) == 1
        return first(config["pn_orders"])
    else
        throw(ArgumentError("Expected PN order `3.5` in the config, or a config with exactly one PN order."))
    end
end

function build_labels()
    return Dict(
        "mc" => L"(1+z)\mathcal{M}",
        "invq" => L"1/q",
        "iota" => L"\iota",
        "chi_eff" => L"\chi_{\mathrm{eff}}",
        "z" => L"z",
    )
end

function xtick_spec(param::String)
    if param == "invq"
        return ([0.0, 0.25, 0.5, 0.75, 1.0], ["0.0", "", "0.5", "", "1.0"])
    elseif param == "iota"
        return ([0.0, π / 4, π / 2, 3π / 4, π], ["0", "", "π/2", "", "π"])
    elseif param == "z"
        return ([0.0, 0.25, 0.5], ["0.00", "0.25", "0.50"])
    else
        return nothing
    end
end

function mc_tick_spec()
    return ([5.0, 20.0, 35.0, 50.0, 65.0, 80.0], ["5", "20", "35", "50", "65", "80"])
end

function build_top_histogram_panel(values::Vector{Float64};
    xlabel, ylabel="", xlim=nothing, xticks_spec=nothing, show_yticks::Bool=true,
    show_ytick_labels::Bool=true, show_xtick_labels::Bool=true,
    left_margin_mm::Real=6, right_margin_mm::Real=4)

    plt = histogram(
        values;
        bins=30,
        normalize=:pdf,
        color=CONTOUR_FILL_COLORS[2],
        linecolor=CONTOUR_LINE_COLOR,
        linewidth=1.0,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
        grid=false,
        framestyle=:box,
        tick_direction=:out,
        dpi=200,
        left_margin=left_margin_mm * Plots.mm,
        right_margin=right_margin_mm * Plots.mm,
        bottom_margin=2Plots.mm,
        top_margin=3Plots.mm,
    )
    xlim === nothing || xlims!(plt, xlim)
    xticks_spec === nothing || xticks!(plt, xticks_spec...)
    if !show_xtick_labels
        if xticks_spec !== nothing
            xticks!(plt, first(xticks_spec), fill("", length(first(xticks_spec))))
        else
            xmin, xmax = Plots.xlims(plt)
            current_ticks = collect(range(xmin, xmax; length=4))
            xticks!(plt, current_ticks, fill("", length(current_ticks)))
        end
    end
    if !show_yticks
        yticks!(plt, (Float64[], String[]))
    elseif !show_ytick_labels
        ymin, ymax = Plots.ylims(plt)
        current_ticks = collect(range(ymin, ymax; length=4))
        yticks!(plt, current_ticks, fill("", length(current_ticks)))
    end
    return plt
end

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#

"""
    read_catalog(fisher_results_file::AbstractString)

Load the event catalog stored in `fisher_results_file` and return it as a
dictionary of event-aligned vectors.

The function reads the datasets stored under `/bbh_catalog/parameter` that are
needed by this plotting script: `"mc"`, `"eta"`, `"chi_1"`, `"chi_2"`,
`"iota"`, and `"z"`.
"""
function read_catalog(fisher_results_file::AbstractString)
    catalog = Dict{String, Vector{Float64}}()

    h5open(fisher_results_file, "r") do file
        catalog_group = file["bbh_catalog"]["parameter"]
        for name in ("mc", "eta", "chi_1", "chi_2", "iota", "z")
            catalog[name] = read(catalog_group, name)
        end
    end

    return catalog
end

"""
    read_results(fisher_results_file::AbstractString, config::Dict, pno::String)

Load the per-waveform Fisher-analysis outputs for post-Newtonian order `pno`
from `fisher_results_file`.

The returned dictionary is keyed by waveform family name (`"PhenomD"` and
`"PhenomHM"`). Each waveform entry contains the event-aligned vectors `"snr"`,
`"isnr"`, `"delta_k"`, and `"invc"` read from the group
`/<pn_order>/<network>/<waveform>`, where the PN-order group name is derived via
`createSED.pnoString(pno)`.
"""
function read_results(fisher_results_file::AbstractString, config::Dict, pno::String)
    results = Dict{String, Dict{String, Vector}}()
    pno_group_name = createSED.pnoString(pno)

    h5open(fisher_results_file, "r") do file
        for wf_name in ("PhenomD", "PhenomHM")
            wf_group = file[pno_group_name][config["network"]][wf_name]
            results[wf_name] = Dict(
                "snr" => read(wf_group, "snr"),
                "isnr" => read(wf_group, "isnr"),
                "delta_k" => read(wf_group, "delta_k"),
                "invc" => read(wf_group, "invc"),
            )
        end
    end

    return results
end

"""
    extend_catalog!(catalog::Dict{String, Vector{Float64}})

Compute derived source quantities from the catalog's chirp mass, symmetric mass
ratio, and aligned spins, and store them back into `catalog`.

The function adds the following event-aligned fields:
- `"m1"`: primary component mass
- `"m2"`: secondary component mass
- `"invq"`: inverse mass ratio `m2 / m1`
- `"chi_eff"`: effective aligned spin

It returns the updated `catalog` dictionary.
"""
function extend_catalog!(catalog::Dict{String, Vector{Float64}})
    mass_ratio_term = sqrt.(max.(1.0 .- 4.0 .* catalog["eta"], 0.0))
    total_mass = catalog["mc"] ./ (catalog["eta"] .^ (3.0 / 5.0))
    m1 = 0.5 .* total_mass .* (1.0 .+ mass_ratio_term)
    m2 = 0.5 .* total_mass .* (1.0 .- mass_ratio_term)

    catalog["m1"] = m1
    catalog["m2"] = m2
    catalog["invq"] = m2 ./ m1
    catalog["chi_eff"] = (m1 .* catalog["chi_1"] .+ m2 .* catalog["chi_2"]) ./ total_mass

    return catalog
end

"""
    build_selection_indices(results::Dict{String, Dict{String, Vector}}, config::Dict)

Construct the event-selection masks used to classify injections from the
per-waveform Fisher-analysis results.

The returned dictionary contains the following event-aligned `BitVector`s:
- `"both_fisher_selected"`: selected by the Fisher/SNR cut for both waveform families
- `"both_observable"`: selected by the full observable cut for both waveform families
- `"potentially_problematic"`: Fisher-selected for both waveforms but not observable for both
- `"d_only"`: observable with `PhenomD` only
- `"hm_only"`: observable with `PhenomHM` only

A waveform-level event is marked as `"fisher_selected"` when its network SNR is
above `config["snr_threshold"]`, the Fisher matrix is invertible, and `delta_k`
is finite and positive. A waveform-level event is marked as `"observable"` when
it is Fisher-selected and its inspiral SNR is above
`config["snr_inspiral_threshold"]`.

The combined masks are then formed across the two waveform families. The
`"d_only"` and `"hm_only"` classes are derived from the waveform-level
`"observable"` masks.
"""
function build_selection_indices(results::Dict{String, Dict{String, Vector}}, config::Dict)
    fisher_valid = Dict{String, BitVector}()
    fisher_selected = Dict{String, BitVector}()
    observable = Dict{String, BitVector}()

    for wf_name in ("PhenomD", "PhenomHM")
        delta_k = Float64.(results[wf_name]["delta_k"])
        snr = Float64.(results[wf_name]["snr"])
        isnr = Float64.(results[wf_name]["isnr"])
        invc = BitVector(results[wf_name]["invc"])

        fisher_valid[wf_name] = BitVector(invc .& isfinite.(delta_k) .& (delta_k .> 0.0))
        fisher_selected[wf_name] = BitVector((snr .> config["snr_threshold"]) .& fisher_valid[wf_name])
        observable[wf_name] = BitVector(fisher_selected[wf_name] .& (isnr .> config["snr_inspiral_threshold"]))
    end

    idx_bfisher = fisher_selected["PhenomD"] .& fisher_selected["PhenomHM"]
    idx_bobs = observable["PhenomD"] .& observable["PhenomHM"]
    idx_prob = idx_bfisher .& .!idx_bobs
    idx_do = observable["PhenomD"] .& .!observable["PhenomHM"]
    idx_hmo = .!observable["PhenomD"] .& observable["PhenomHM"]

    return Dict(
        "PhenomD_fisher_selected" => fisher_selected["PhenomD"],
        "PhenomHM_fisher_selected" => fisher_selected["PhenomHM"],
        "PhenomD_observable" => observable["PhenomD"],
        "PhenomHM_observable" => observable["PhenomHM"],
        "both_fisher_selected" => idx_bfisher,
        "both_observable" => idx_bobs,
        "potentially_problematic" => idx_prob,
        "d_only" => idx_do,
        "hm_only" => idx_hmo,
    )
end

function delta_ratio(results::Dict{String, Dict{String, Vector}})
    delta_hm = Float64.(results["PhenomHM"]["delta_k"])
    delta_d = Float64.(results["PhenomD"]["delta_k"])
    return log10.(delta_hm ./ delta_d)
end

"""
    contour_levels_from_density(density::AbstractMatrix{<:Real},
                                enclosed_masses=(0.98, 0.90, 0.30))

Compute contour thresholds for a two-dimensional density estimate so that each
threshold corresponds to a requested enclosed probability mass.

The density values are flattened, sorted from highest to lowest, and cumulatively
integrated. For each entry in `enclosed_masses`, the returned threshold is the
density value at which the cumulative enclosed mass first reaches that target.
An additional upper level slightly above the maximum density is appended so the
result can be passed directly to filled-contour plotting routines.
"""
function contour_levels_from_density(density::AbstractMatrix{<:Real}, enclosed_masses=(0.98, 0.90, 0.30))
    weights = vec(Float64.(density))
    total = sum(weights)
    total == 0.0 && return [0.0, 1.0]

    sorted_weights = sort(weights; rev=true)
    cumulative = cumsum(sorted_weights) ./ total

    function level_for_mass(target_mass)
        idx = searchsortedfirst(cumulative, target_mass)
        idx = clamp(idx, 1, length(sorted_weights))
        return sorted_weights[idx]
    end

    levels = [max(0.0, level_for_mass(mass)) for mass in enclosed_masses]
    return vcat(levels, maximum(sorted_weights) + eps())
end

"""
    kde_catalog_contours(catalog::Dict{String, Vector{Float64}}, x_key::String;
                         z_threshold::Float64=0.5, npoints::Int=70,
                         enclosed_masses=(0.98, 0.90, 0.30))

Perform a two-dimensional kernel-density estimate for the catalog distribution
of `catalog[x_key]` versus chirp mass `catalog["mc"]`, using only events with
redshift smaller than `z_threshold`.

Returns `(x_grid, y_grid, density, levels)`, where `x_grid` and `y_grid` are
the KDE evaluation grids, `density` is the estimated probability density on the
grid, and `levels` are contour thresholds corresponding to the requested
enclosed probability masses.
"""
function kde_catalog_contours(catalog::Dict{String, Vector{Float64}}, x_key::String;
    z_threshold::Float64=FIG9_Z_THRESHOLD, npoints::Int=70, enclosed_masses=(0.98, 0.90, 0.30))

    haskey(catalog, x_key) || throw(ArgumentError("Catalog does not contain key: $(x_key)"))

    selected = BitVector(catalog["z"] .< z_threshold)
    x = Float64.(catalog[x_key][selected])
    y = Float64.(catalog["mc"][selected])

    isempty(x) && throw(ArgumentError("No catalog events satisfy z < $(z_threshold)."))

    kde_result = kde((x, y); npoints=(npoints, npoints))
    density = Float64.(kde_result.density)
    levels = contour_levels_from_density(density, enclosed_masses)

    return kde_result.x, kde_result.y, density, levels
end

function build_contour_panel(catalog::Dict{String, Vector{Float64}}, x_key::String;
    xlabel, ylabel="", xlim=nothing, show_yticks::Bool=true, show_ytick_labels::Bool=true,
    left_margin_mm::Real=6, right_margin_mm::Real=4, z_threshold::Float64=FIG9_Z_THRESHOLD)

    x_grid, y_grid, density, levels = kde_catalog_contours(catalog, x_key; z_threshold=z_threshold)

    plt = contourf(
        x_grid,
        y_grid,
        density';
        levels=levels,
        c=CONTOUR_FILL_COLORS,
        linecolor=false,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
        grid=false,
        framestyle=:box,
        dpi=200,
        left_margin=left_margin_mm * Plots.mm,
        right_margin=right_margin_mm * Plots.mm,
        bottom_margin=10Plots.mm,
        top_margin=1Plots.mm,
    )
    contour!(
        plt,
        x_grid,
        y_grid,
        density';
        levels=levels[1:end-1],
        color=CONTOUR_LINE_COLOR,
        linewidth=1.2,
        alpha=0.8,
        label=false,
    )
    xlim === nothing || xlims!(plt, xlim)
    ylims!(plt, FIG9_MC_LIMITS)
    tick_spec = xtick_spec(x_key)
    tick_spec === nothing || xticks!(plt, tick_spec...)
    if !show_yticks
        yticks!(plt, (Float64[], String[]))
    else
        y_tick_spec = mc_tick_spec()
        if show_ytick_labels
            yticks!(plt, y_tick_spec...)
        else
            yticks!(plt, first(y_tick_spec), fill("", length(first(y_tick_spec))))
        end
    end
    return plt
end

function build_colorbar_panel(color_lims)
    midpoint = 0.5 * (color_lims[1] + color_lims[2])
    plt = scatter(
        [0.0],
        [0.0];
        marker_z=[midpoint],
        color=OBSERVABLE_COLORMAP,
        clims=color_lims,
        markersize=0,
        markerstrokewidth=0,
        alpha=0.0,
        colorbar=true,
        colorbar_title=IMPROVEMENT_COLORBAR_TITLE,
        framestyle=:none,
        grid=false,
        xticks=false,
        yticks=false,
        xlims=(0.0, 1.0),
        ylims=(0.0, 1.0),
        dpi=200,
        left_margin=4Plots.mm,
        right_margin=12Plots.mm,
        bottom_margin=10Plots.mm,
        top_margin=1Plots.mm,
        label=false,
        foreground_color_subplot=:white,
        background_color_subplot=:white,
    )
    return plt
end

function problematic_box_halfwidth(xlim::Tuple{<:Real, <:Real}; halfheight::Float64=1.5)
    chirp_mass_range = 80.0 - 5.0
    panel_height_over_width = 1.6
    return halfheight * panel_height_over_width * (Float64(xlim[2]) - Float64(xlim[1])) / chirp_mass_range
end

function add_problematic_boxes!(plt, x::Vector{Float64}, y::Vector{Float64}, xlim::Tuple{<:Real, <:Real})
    halfheight = 0.9
    halfwidth = problematic_box_halfwidth(xlim; halfheight=halfheight)

    for idx in eachindex(x)
        box = Shape(
            [x[idx] - halfwidth, x[idx] + halfwidth, x[idx] + halfwidth, x[idx] - halfwidth],
            [y[idx] - halfheight, y[idx] - halfheight, y[idx] + halfheight, y[idx] + halfheight],
        )
        plot!(
            plt,
            box;
            fillalpha=0.0,
            linecolor=:black,
            linewidth=1.4,
            label=false,
        )
    end

    return plt
end

function add_fisher_selected_scatter!(plt, x::Vector{Float64}, y::Vector{Float64}, ratio::Vector{Float64},
    selected::BitVector, fisher_selected::BitVector, problematic::BitVector, color_lims,
    xlim::Tuple{<:Real, <:Real}; show_colorbar::Bool=false)

    idx = BitVector(selected .& fisher_selected)
    idx_problematic = BitVector(selected .& problematic)
    scatter!(
        plt,
        x[idx],
        y[idx];
        marker_z=ratio[idx],
        color=OBSERVABLE_COLORMAP,
        clims=color_lims,
        markersize=6.75,
        markerstrokewidth=0.0,
        alpha=0.95,
        colorbar=show_colorbar,
        colorbar_title=show_colorbar ? IMPROVEMENT_COLORBAR_TITLE : "",
        label=false,
    )
    add_problematic_boxes!(plt, x[idx_problematic], y[idx_problematic], xlim)

    return plt
end

function build_plot(catalog::Dict{String, Vector{Float64}}, results::Dict{String, Dict{String, Vector}}, indices::Dict{String, BitVector})
    extend_catalog!(catalog)
    labels = build_labels()
    kde_selected = BitVector(catalog["z"] .< FIG9_Z_THRESHOLD)
    selected = BitVector(kde_selected .& indices["both_fisher_selected"])
    ratio = delta_ratio(results)
    valid_ratio = ratio[selected]
    color_lims = isempty(valid_ratio) ? (-1.0, 1.0) : (minimum(valid_ratio), maximum(valid_ratio))

    h1 = build_top_histogram_panel(catalog["invq"][kde_selected]; xlabel="", ylabel="", xlim=(0.0, 1.0), xticks_spec=xtick_spec("invq"), show_yticks=false, show_xtick_labels=false, left_margin_mm=12, right_margin_mm=5)
    h2 = build_top_histogram_panel(catalog["iota"][kde_selected]; xlabel="", ylabel="", xlim=(0.0, π), xticks_spec=xtick_spec("iota"), show_yticks=false, show_xtick_labels=false, left_margin_mm=4, right_margin_mm=4)
    h3 = build_top_histogram_panel(catalog["chi_eff"][kde_selected]; xlabel="", ylabel="", xlim=(-1.0, 1.0), show_yticks=false, show_xtick_labels=false, left_margin_mm=4, right_margin_mm=4)
    h4 = build_top_histogram_panel(catalog["z"][kde_selected]; xlabel="", ylabel="", xlim=(0.0, FIG9_Z_THRESHOLD), xticks_spec=xtick_spec("z"), show_yticks=false, show_xtick_labels=false, left_margin_mm=4, right_margin_mm=4)
    h5 = build_top_histogram_panel(catalog["mc"][kde_selected]; xlabel=labels["mc"], ylabel="", xlim=(5.0, 80.0), xticks_spec=mc_tick_spec(), show_yticks=false, left_margin_mm=4, right_margin_mm=10)

    p1 = build_contour_panel(catalog, "invq"; xlabel=labels["invq"], ylabel=labels["mc"], xlim=(0.0, 1.0), show_yticks=true, left_margin_mm=12, right_margin_mm=5, z_threshold=FIG9_Z_THRESHOLD)
    add_fisher_selected_scatter!(p1, catalog["invq"], catalog["mc"], ratio, selected, indices["both_fisher_selected"], indices["potentially_problematic"], color_lims, (0.0, 1.0))
    ylims!(p1, FIG9_MC_LIMITS)
    p2 = build_contour_panel(catalog, "iota"; xlabel=labels["iota"], xlim=(0.0, π), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=4, z_threshold=FIG9_Z_THRESHOLD)
    add_fisher_selected_scatter!(p2, catalog["iota"], catalog["mc"], ratio, selected, indices["both_fisher_selected"], indices["potentially_problematic"], color_lims, (0.0, π))
    ylims!(p2, FIG9_MC_LIMITS)
    p3 = build_contour_panel(catalog, "chi_eff"; xlabel=labels["chi_eff"], xlim=(-1.0, 1.0), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=4, z_threshold=FIG9_Z_THRESHOLD)
    add_fisher_selected_scatter!(p3, catalog["chi_eff"], catalog["mc"], ratio, selected, indices["both_fisher_selected"], indices["potentially_problematic"], color_lims, (-1.0, 1.0))
    ylims!(p3, FIG9_MC_LIMITS)
    p4 = build_contour_panel(catalog, "z"; xlabel=labels["z"], xlim=(0.0, FIG9_Z_THRESHOLD), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=4, z_threshold=FIG9_Z_THRESHOLD)
    add_fisher_selected_scatter!(p4, catalog["z"], catalog["mc"], ratio, selected, indices["both_fisher_selected"], indices["potentially_problematic"], color_lims, (0.0, FIG9_Z_THRESHOLD))
    ylims!(p4, FIG9_MC_LIMITS)
    p5 = build_colorbar_panel(color_lims)

    return plot(
        h1, h2, h3, h4, h5,
        p1, p2, p3, p4, p5;
        layout=grid(2, 5, heights=[0.34, 0.66]),
        size=(1900, 900),
        left_margin=6Plots.mm,
        right_margin=6Plots.mm,
        bottom_margin=8Plots.mm,
        top_margin=2Plots.mm,
    )
end

function run_pop_prob(config::Dict)
    fisher_results_file = get_fisher_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    pno = get_target_pno(config)
    catalog = read_catalog(fisher_results_file)
    results = read_results(fisher_results_file, config, pno)
    indices = build_selection_indices(results, config)
    observable_redshifts = catalog["z"][indices["both_observable"]]
    fisher_selected_redshifts = catalog["z"][indices["both_fisher_selected"]]
    max_observable_z = isempty(observable_redshifts) ? NaN : maximum(observable_redshifts)
    max_fisher_selected_z = isempty(fisher_selected_redshifts) ? NaN : maximum(fisher_selected_redshifts)
    println("Maximum redshift among observable events: $(max_observable_z)")
    println("Maximum redshift among Fisher-selected events: $(max_fisher_selected_z)")

    plt = build_plot(catalog, results, indices)

    png_output_file, pdf_output_file = get_plot_output_files(config)
    mkpath(dirname(png_output_file))
    savefig(plt, png_output_file)
    savefig(plt, pdf_output_file)
    sleep(5)

    println("Saved figure to $(png_output_file)")
    println("Saved figure to $(pdf_output_file)")
    return png_output_file
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
function main(args=ARGS)
    if length(args) != 1
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml>"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    config = read_config(config_file)

    return run_pop_prob(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
