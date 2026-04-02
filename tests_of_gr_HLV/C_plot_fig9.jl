using TOML
using HDF5
using KernelDensity
using CairoMakie
using LaTeXStrings
using Colors

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const CONTOUR_FILL_COLORS = [:aliceblue, :lightblue, :cornflowerblue]
const CONTOUR_LINE_COLOR = :royalblue4
const OBSERVABLE_COLORMAP = :viridis
const FIG9_Z_THRESHOLD = 0.5
const FIG9_MC_LIMITS = (5.0, 80.0)
const FIG9_SIZE = (1800, 760)
const FIG9_HIST_FRACTION = 0.32
const FIG9_COL_GAP = 6*5
const FIG9_ROW_GAP = 4*5
const FIG9_LABELSIZE =28
const HIST_BASE_COLOR = :lightblue
const HIST_FISHER_COLOR = :crimson
const HIST_OBSERVABLE_COLOR = :darkorange

function fig9_side_width_fraction()
    return FIG9_HIST_FRACTION * FIG9_SIZE[2] / FIG9_SIZE[1]
end

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

function kde_catalog_contours(catalog::Dict{String, Vector{Float64}}, x_key::String;
    z_threshold::Float64=FIG9_Z_THRESHOLD, npoints::Int=70, enclosed_masses=(0.98, 0.90, 0.30))

    selected = BitVector(catalog["z"] .< z_threshold)
    x = Float64.(catalog[x_key][selected])
    y = Float64.(catalog["mc"][selected])
    isempty(x) && throw(ArgumentError("No catalog events satisfy z < $(z_threshold)."))

    kde_result = kde((x, y); npoints=(npoints, npoints))
    density = Float64.(kde_result.density)
    levels = contour_levels_from_density(density, enclosed_masses)

    return kde_result.x, kde_result.y, density, levels
end

function histogram_bin_counts(values::Vector{Float64}, limits::Tuple{<:Real, <:Real}; nbins::Int=20)
    xmin = Float64(limits[1])
    xmax = Float64(limits[2])
    edges = collect(range(xmin, xmax; length=nbins + 1))
    counts = zeros(Float64, nbins)

    for value in values
        if !isfinite(value) || value < xmin || value > xmax
            continue
        end
        idx = value == xmax ? nbins : searchsortedlast(edges, value)
        idx = clamp(idx, 1, nbins)
        counts[idx] += 1.0
    end

    return edges, counts
end

function histogram_bin_densities(values::Vector{Float64}, limits::Tuple{<:Real, <:Real}; nbins::Int=20)
    edges, counts = histogram_bin_counts(values, limits; nbins=nbins)
    bin_width = (Float64(limits[2]) - Float64(limits[1])) / nbins
    total = sum(counts)
    if total > 0 && bin_width > 0
        counts ./= (total * bin_width)
    end
    return edges, counts
end

function scale_density_to_fisher_peak(density_values::Vector{Float64}, fisher_peak::Float64; factor::Float64=1.3)
    density_peak = maximum(density_values)
    target_peak = factor * fisher_peak
    if density_peak > 0 && target_peak > 0
        return density_values .* (target_peak / density_peak)
    end
    return zeros(length(density_values))
end

function x_ticks_for_param(param::String)
    if param == "invq"
        return ([0.0, 0.5, 1.0], ["0.0", "0.5", "1.0"])
    elseif param == "iota"
        return ([0.0, π / 2, π], ["0", "π/2", "π"])
    elseif param == "z"
        return ([0.0, 0.25, 0.5], ["0.00", "0.25", "0.50"])
    else
        return nothing
    end
end

function mc_ticks()
    return ([5.0, 20.0, 35.0, 50.0, 65.0, 80.0], ["5", "20", "35", "50", "65", "80"])
end

function hide_hist_axis!(ax::Axis; bottom_spine::Bool=true, left_spine::Bool=false, hide_x::Bool=true, hide_y::Bool=true)
    hidedecorations!(ax; label=hide_x && hide_y, ticklabels=hide_x && hide_y, ticks=hide_x && hide_y)
    if hide_x
        hidexdecorations!(ax; label=true, ticklabels=true, ticks=true)
    end
    if hide_y
        hideydecorations!(ax; label=true, ticklabels=true, ticks=true)
    end
    if bottom_spine
        hidespines!(ax, :l, :r, :t)
    elseif left_spine
        hidespines!(ax, :r, :t, :b)
    else
        hidespines!(ax)
    end
    return ax
end

function step_xy(edges::Vector{Float64}, counts::Vector{Float64})
    xcoords = Float64[]
    ycoords = Float64[]
    push!(xcoords, edges[1]); push!(ycoords, 0.0)
    for i in eachindex(counts)
        push!(xcoords, edges[i]); push!(ycoords, counts[i])
        push!(xcoords, edges[i + 1]); push!(ycoords, counts[i])
    end
    push!(xcoords, edges[end]); push!(ycoords, 0.0)
    return xcoords, ycoords
end

function side_step_xy(edges::Vector{Float64}, counts::Vector{Float64})
    xcoords = Float64[]
    ycoords = Float64[]
    push!(xcoords, 0.0); push!(ycoords, edges[1])
    for i in eachindex(counts)
        push!(xcoords, counts[i]); push!(ycoords, edges[i])
        push!(xcoords, counts[i]); push!(ycoords, edges[i + 1])
    end
    push!(xcoords, 0.0); push!(ycoords, edges[end])
    return xcoords, ycoords
end

function draw_top_histogram!(ax::Axis, base_values::Vector{Float64}, fisher_values::Vector{Float64},
    observable_values::Vector{Float64}, xlim::Tuple{<:Real, <:Real}; nbins::Int=20)
    edges, base_density = histogram_bin_densities(base_values, xlim; nbins=nbins)
    _, fisher_counts = histogram_bin_counts(fisher_values, xlim; nbins=nbins)
    _, observable_counts = histogram_bin_counts(observable_values, xlim; nbins=nbins)
    base_counts = scale_density_to_fisher_peak(base_density, maximum(fisher_counts))
    ymax = maximum(vcat(base_counts, fisher_counts, observable_counts))
    ymax = ymax > 0 ? 1.05 * ymax : 1.0
    xlims!(ax, Float64(xlim[1]), Float64(xlim[2]))
    ylims!(ax, 0.0, ymax)

    for i in eachindex(base_counts)
        points = Point2f[
            (edges[i], 0.0),
            (edges[i + 1], 0.0),
            (edges[i + 1], base_counts[i]),
            (edges[i], base_counts[i]),
        ]
        poly!(ax, points; color=HIST_BASE_COLOR, strokecolor=:transparent)
    end

    fisher_x, fisher_y = step_xy(edges, fisher_counts)
    observable_x, observable_y = step_xy(edges, observable_counts)
    lines!(ax, fisher_x, fisher_y; color=HIST_FISHER_COLOR, linewidth=2)
    lines!(ax, observable_x, observable_y; color=HIST_OBSERVABLE_COLOR, linewidth=2)
    lines!(ax, [Float64(xlim[1]), Float64(xlim[2])], [0.0, 0.0]; color=:black, linewidth=1)
    return ax
end

function draw_side_histogram!(ax::Axis, base_values::Vector{Float64}, fisher_values::Vector{Float64},
    observable_values::Vector{Float64}; nbins::Int=20)
    edges, base_density = histogram_bin_densities(base_values, FIG9_MC_LIMITS; nbins=nbins)
    _, fisher_counts = histogram_bin_counts(fisher_values, FIG9_MC_LIMITS; nbins=nbins)
    _, observable_counts = histogram_bin_counts(observable_values, FIG9_MC_LIMITS; nbins=nbins)
    base_counts = scale_density_to_fisher_peak(base_density, maximum(fisher_counts))
    xmax = maximum(vcat(base_counts, fisher_counts, observable_counts))
    xmax = xmax > 0 ? 1.05 * xmax : 1.0
    xlims!(ax, 0.0, xmax)
    ylims!(ax, FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2])

    for i in eachindex(base_counts)
        points = Point2f[
            (0.0, edges[i]),
            (base_counts[i], edges[i]),
            (base_counts[i], edges[i + 1]),
            (0.0, edges[i + 1]),
        ]
        poly!(ax, points; color=HIST_BASE_COLOR, strokecolor=:transparent)
    end

    fisher_x, fisher_y = side_step_xy(edges, fisher_counts)
    observable_x, observable_y = side_step_xy(edges, observable_counts)
    lines!(ax, fisher_x, fisher_y; color=HIST_FISHER_COLOR, linewidth=2)
    lines!(ax, observable_x, observable_y; color=HIST_OBSERVABLE_COLOR, linewidth=2)
    lines!(ax, [0.0, 0.0], [FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2]]; color=:black, linewidth=1)
    return ax
end

function problematic_box_halfwidth(xlim::Tuple{<:Real, <:Real}; halfheight::Float64=1.5)
    chirp_mass_range = FIG9_MC_LIMITS[2] - FIG9_MC_LIMITS[1]
    panel_height_over_width = 1.6
    return halfheight * panel_height_over_width * (Float64(xlim[2]) - Float64(xlim[1])) / chirp_mass_range
end

function add_problematic_boxes!(ax::Axis, x::Vector{Float64}, y::Vector{Float64}, xlim::Tuple{<:Real, <:Real})
    halfheight = 0.9
    halfwidth = problematic_box_halfwidth(xlim; halfheight=halfheight)

    for idx in eachindex(x)
        points = Point2f[
            (x[idx] - halfwidth, y[idx] - halfheight),
            (x[idx] + halfwidth, y[idx] - halfheight),
            (x[idx] + halfwidth, y[idx] + halfheight),
            (x[idx] - halfwidth, y[idx] + halfheight),
        ]
        poly!(ax, points; color=:transparent, strokecolor=:black, strokewidth=1.4)
    end
    return ax
end

function build_main_panel!(ax::Axis, catalog::Dict{String, Vector{Float64}}, results::Dict{String, Dict{String, Vector}},
    indices::Dict{String, BitVector}, x_key::String, xlim::Tuple{<:Real, <:Real}, ratio::Vector{Float64}, color_lims;
    xlabel, ylabel="", show_yticks::Bool=true, show_yticklabels::Bool=true)

    x_grid, y_grid, density, levels = kde_catalog_contours(catalog, x_key; z_threshold=FIG9_Z_THRESHOLD)

    contourf!(ax, x_grid, y_grid, density; levels=levels, colormap=CONTOUR_FILL_COLORS)
    contour!(ax, x_grid, y_grid, density; levels=levels[1:end-1], color=CONTOUR_LINE_COLOR, linewidth=1.2)

    selected = BitVector((catalog["z"] .< FIG9_Z_THRESHOLD) .& indices["both_fisher_selected"])
    problematic = BitVector((catalog["z"] .< FIG9_Z_THRESHOLD) .& indices["potentially_problematic"])

    xvals = catalog[x_key][selected]
    yvals = catalog["mc"][selected]
    cvals = ratio[selected]
    scatter!(ax, xvals, yvals; color=cvals, colormap=OBSERVABLE_COLORMAP, colorrange=color_lims, markersize=12)
    add_problematic_boxes!(ax, catalog[x_key][problematic], catalog["mc"][problematic], xlim)

    xlims!(ax, Float64(xlim[1]), Float64(xlim[2]))
    ylims!(ax, FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2])
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.xlabelsize = FIG9_LABELSIZE
    ax.ylabelsize = FIG9_LABELSIZE
    ax.xtickalign = 1
    ax.ytickalign = 1

    tick_spec = x_ticks_for_param(x_key)
    if tick_spec !== nothing
        ax.xticks = tick_spec
    end

    if show_yticks
        ax.yticks = mc_ticks()
    else
        hideydecorations!(ax; ticks=true, ticklabels=true, label=false)
    end
    if show_yticks && !show_yticklabels
        hideydecorations!(ax; ticklabels=true, ticks=false, label=false)
    end

    return ax
end

function build_plot(catalog::Dict{String, Vector{Float64}}, results::Dict{String, Dict{String, Vector}}, indices::Dict{String, BitVector})
    CairoMakie.activate!()
    extend_catalog!(catalog)
    labels = build_labels()
    kde_selected = BitVector(catalog["z"] .< FIG9_Z_THRESHOLD)
    ratio = delta_ratio(results)
    fisher_hist_selected = BitVector(kde_selected .& indices["both_fisher_selected"])
    observable_hist_selected = BitVector(kde_selected .& indices["both_observable"])
    scatter_selected = BitVector(kde_selected .& indices["both_fisher_selected"])
    valid_ratio = ratio[scatter_selected]
    color_lims = isempty(valid_ratio) ? (-1.0, 1.0) : (minimum(valid_ratio), maximum(valid_ratio))

    fig = Figure(size=FIG9_SIZE, backgroundcolor=:white)
    main_width = (1.0 - fig9_side_width_fraction()) / 4.0

    top_axes = [Axis(fig[1, i], backgroundcolor=:transparent) for i in 1:4]
    blank_ax = Axis(fig[1, 5], backgroundcolor=:transparent)
    hide_hist_axis!(blank_ax; bottom_spine=false)

    main_axes = [Axis(fig[2, i], backgroundcolor=:white) for i in 1:4]
    side_ax = Axis(fig[2, 5], backgroundcolor=:white)

    for col in 1:4
        colsize!(fig.layout, col, Relative(main_width))
    end
    colsize!(fig.layout, 5, Relative(fig9_side_width_fraction()))
    rowsize!(fig.layout, 1, Relative(FIG9_HIST_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG9_HIST_FRACTION))
    colgap!(fig.layout, FIG9_COL_GAP)
    rowgap!(fig.layout, FIG9_ROW_GAP)

    draw_top_histogram!(top_axes[1], catalog["invq"][kde_selected], catalog["invq"][fisher_hist_selected], catalog["invq"][observable_hist_selected], (0.0, 1.0))
    draw_top_histogram!(top_axes[2], catalog["iota"][kde_selected], catalog["iota"][fisher_hist_selected], catalog["iota"][observable_hist_selected], (0.0, π))
    draw_top_histogram!(top_axes[3], catalog["chi_eff"][kde_selected], catalog["chi_eff"][fisher_hist_selected], catalog["chi_eff"][observable_hist_selected], (-1.0, 1.0))
    draw_top_histogram!(top_axes[4], catalog["z"][kde_selected], catalog["z"][fisher_hist_selected], catalog["z"][observable_hist_selected], (0.0, FIG9_Z_THRESHOLD))
    foreach(ax -> hide_hist_axis!(ax; bottom_spine=true), top_axes)

    build_main_panel!(main_axes[1], catalog, results, indices, "invq", (0.0, 1.0), ratio, color_lims;
        xlabel=labels["invq"], ylabel=labels["mc"], show_yticks=true, show_yticklabels=true)
    build_main_panel!(main_axes[2], catalog, results, indices, "iota", (0.0, π), ratio, color_lims;
        xlabel=labels["iota"], show_yticks=true, show_yticklabels=false)
    build_main_panel!(main_axes[3], catalog, results, indices, "chi_eff", (-1.0, 1.0), ratio, color_lims;
        xlabel=labels["chi_eff"], show_yticks=true, show_yticklabels=false)
    build_main_panel!(main_axes[4], catalog, results, indices, "z", (0.0, FIG9_Z_THRESHOLD), ratio, color_lims;
        xlabel=labels["z"], show_yticks=true, show_yticklabels=false)
    linkyaxes!(main_axes...)

    draw_side_histogram!(side_ax, catalog["mc"][kde_selected], catalog["mc"][fisher_hist_selected], catalog["mc"][observable_hist_selected])
    hide_hist_axis!(side_ax; left_spine=true, hide_x=true, hide_y=true)
    ylims!(side_ax, FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2])

    return fig
end

function run_plot_fig9_makie(config::Dict)
    fisher_results_file = get_fisher_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    pno = get_target_pno(config)
    catalog = read_catalog(fisher_results_file)
    results = read_results(fisher_results_file, config, pno)
    indices = build_selection_indices(results, config)

    fig = build_plot(catalog, results, indices)

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

    return run_plot_fig9_makie(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
