using TOML
using HDF5
using KernelDensity
using Plots
using LaTeXStrings

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const CONTOUR_FILL_COLORS = [:aliceblue, :lightblue, :cornflowerblue]
const CONTOUR_LINE_COLOR = :royalblue4

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_plot_output_files(config::Dict)
    output_dir = joinpath(@__DIR__, config["plot_outdir"])
    stem = "pop_prob_$(config["plot_tag"])"
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
- `"both_sufficient_snr"`: above the network SNR threshold for both waveform families
- `"both_observable"`: selected by the full observable cut for both waveform families
- `"potentially_problematic"`: sufficient SNR for both waveforms but not observable for both
- `"d_only"`: observable with `PhenomD` only
- `"hm_only"`: observable with `PhenomHM` only

A waveform-level event is marked as `"sufficient_snr"` when its network SNR is
above `config["snr_threshold"]`. A waveform-level event is marked as
`"observable"` when the Fisher matrix is invertible, `delta_k` is finite and
positive, and its inspiral SNR is above `config["snr_inspiral_threshold"]`.

The combined masks are then formed across the two waveform families. The
`"d_only"` and `"hm_only"` classes are derived from the waveform-level
`"observable"` masks.
"""
function build_selection_indices(results::Dict{String, Dict{String, Vector}}, config::Dict)
    fisher_valid = Dict{String, BitVector}()
    sufficient_snr = Dict{String, BitVector}()
    observable = Dict{String, BitVector}()

    for wf_name in ("PhenomD", "PhenomHM")
        delta_k = Float64.(results[wf_name]["delta_k"])
        snr = Float64.(results[wf_name]["snr"])
        isnr = Float64.(results[wf_name]["isnr"])
        invc = BitVector(results[wf_name]["invc"])

        fisher_valid[wf_name] = BitVector(invc .& isfinite.(delta_k) .& (delta_k .> 0.0))
        sufficient_snr[wf_name] = BitVector(snr .> config["snr_threshold"])
        observable[wf_name] = BitVector(
            fisher_valid[wf_name] .&
            (isnr .> config["snr_inspiral_threshold"]) .&
            sufficient_snr[wf_name]
        )
    end

    idx_bsnr = sufficient_snr["PhenomD"] .& sufficient_snr["PhenomHM"]
    idx_bobs = observable["PhenomD"] .& observable["PhenomHM"]
    idx_prob = idx_bsnr .& .!idx_bobs
    idx_do = observable["PhenomD"] .& .!observable["PhenomHM"]
    idx_hmo = .!observable["PhenomD"] .& observable["PhenomHM"]

    return Dict(
        "PhenomD_sufficient_snr" => sufficient_snr["PhenomD"],
        "PhenomHM_sufficient_snr" => sufficient_snr["PhenomHM"],
        "PhenomD_observable" => observable["PhenomD"],
        "PhenomHM_observable" => observable["PhenomHM"],
        "both_sufficient_snr" => idx_bsnr,
        "both_observable" => idx_bobs,
        "potentially_problematic" => idx_prob,
        "d_only" => idx_do,
        "hm_only" => idx_hmo,
    )
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
    z_threshold::Float64=0.5, npoints::Int=70, enclosed_masses=(0.98, 0.90, 0.30))

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
    left_margin_mm::Real=6, right_margin_mm::Real=4)

    x_grid, y_grid, density, levels = kde_catalog_contours(catalog, x_key)

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
        top_margin=4Plots.mm,
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

function build_plot(catalog::Dict{String, Vector{Float64}})
    extend_catalog!(catalog)
    labels = build_labels()
    p1 = build_contour_panel(catalog, "invq"; xlabel=labels["invq"], ylabel=labels["mc"], xlim=(0.0, 1.0), show_yticks=true, left_margin_mm=12, right_margin_mm=5)
    p2 = build_contour_panel(catalog, "iota"; xlabel=labels["iota"], xlim=(0.0, π), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=4)
    p3 = build_contour_panel(catalog, "chi_eff"; xlabel=labels["chi_eff"], xlim=(-1.0, 1.0), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=4)
    p4 = build_contour_panel(catalog, "z"; xlabel=labels["z"], xlim=(0.0, 0.5), show_yticks=true, show_ytick_labels=false, left_margin_mm=4, right_margin_mm=10)

    return plot(
        p1, p2, p3, p4;
        layout=(1, 4),
        size=(1600, 650),
        link=:y,
        left_margin=6Plots.mm,
        right_margin=6Plots.mm,
        bottom_margin=8Plots.mm,
        top_margin=4Plots.mm,
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
    sufficient_snr_redshifts = catalog["z"][indices["both_sufficient_snr"]]
    max_observable_z = isempty(observable_redshifts) ? NaN : maximum(observable_redshifts)
    max_sufficient_snr_z = isempty(sufficient_snr_redshifts) ? NaN : maximum(sufficient_snr_redshifts)
    println("Maximum redshift among observable events: $(max_observable_z)")
    println("Maximum redshift among sufficient-SNR events: $(max_sufficient_snr_z)")

    plt = build_plot(catalog)

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
