using TOML
using HDF5
using KernelDensity
using CairoMakie
using LaTeXStrings
using Colors

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

const FIG9_BLUE = "#0072B2"
const FIG9_BLUE_FILL = "#A8CFE5"
const FIG9_ORANGE = "#D55E00"
const FIG9_ORANGE_FILL = "#F1C8A8"
const FIG9_GREEN = "#009E73"
const FIG9_GREEN_FILL = "#A8DCCF"
const FIG9_PURPLE = "#CC79A7"
const CONTOUR_FILL_COLORS = ["#F7F7F7", "#E0E0E0", "#BDBDBD"]
const CONTOUR_LINE_COLOR = "#8F8F8F"
const IMPROVEMENT_COLORMAP = [
    "#E86A17",
    "#D46B24",
    "#BF6B32",
    "#A76943",
    "#8F6556",
    "#75637A",
    "#5F6B98",
    "#4973B6",
    "#347BC7",
    "#1E6CCB",
]
const IMPROVEMENT_COLORS = RGB{Float64}.(parse.(Ref(Colorant), IMPROVEMENT_COLORMAP))
const IMPROVEMENT_COLORBAR_TITLE = L"\log_{10}\!\left(\Delta k_{\mathrm{HM}} / \Delta k_{\mathrm{D}}\right)"
const FIG9_Z_THRESHOLD = 0.5
const FIG9_MC_LIMITS = (5.0, 80.0)
const FIG9_SIZE = (1800, 760)
const FIG9_HIST_FRACTION = 0.32
const FIG9_TOP_HIST_FRACTION = 1.25 * FIG9_HIST_FRACTION
const FIG9_SIDE_HIST_WIDTH_SCALE = 1.40
const FIG9_COL_GAP = 10*5
const FIG9_ROW_GAP = 2*5
const FIG9_AXIS_LABELSIZE = 34
const FIG9_AXIS_TICK_LABELSIZE = 24
const FIG9_COLORBAR_TICK_LABELSIZE = 24
const FIG9_MARKER_SIZE = 18
const FIG9_MARKER_STROKE_WIDTH = 2.2
const FIG9_LEGEND_MARKER_SIZE = 15
const FIG9_LEGEND_LABELSIZE = 22
const FIG9_COLORBAR_LABELSIZE = FIG9_LEGEND_LABELSIZE
const FIG9_LEGEND_COLORBAR_PAD = 14
const HIST_BASE_COLOR = "#D9D9D9"
const HIST_BASE_EDGE_COLOR = "#8F8F8F"
const HIST_LINE_WIDTH = 4
const HIST_OBSERVABLE_LINE_WIDTH = 2
const HIST_OUTLINE_COLOR = :black
const LEGEND_IMPROVEMENT_VALUES = [-2.0, -1.0, 0.0]

function fig9_side_width_fraction()
    return FIG9_SIDE_HIST_WIDTH_SCALE * FIG9_HIST_FRACTION * FIG9_SIZE[2] / FIG9_SIZE[1]
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

function bin_index(value::Real, edges::Vector{Float64}, limits::Tuple{<:Real, <:Real})
    xmin = Float64(limits[1])
    xmax = Float64(limits[2])
    if !isfinite(value) || value < xmin || value > xmax
        return nothing
    end
    nbins = length(edges) - 1
    idx = value == xmax ? nbins : searchsortedlast(edges, Float64(value))
    return clamp(idx, 1, nbins)
end

function histogram_bin_ratios(values::Vector{Float64}, ratios::Vector{Float64}, edges::Vector{Float64},
    limits::Tuple{<:Real, <:Real})

    binned_ratios = [Float64[] for _ in 1:(length(edges) - 1)]
    for (value, ratio) in zip(values, ratios)
        !isfinite(ratio) && continue
        idx = bin_index(value, edges, limits)
        idx === nothing && continue
        push!(binned_ratios[idx], Float64(ratio))
    end
    foreach(bin -> sort!(bin; by=abs), binned_ratios)
    return binned_ratios
end

function ratio_color(value::Real, color_lims::Tuple{<:Real, <:Real})
    lo = Float64(color_lims[1])
    hi = Float64(color_lims[2])
    t = hi == lo ? 0.5 : clamp((Float64(value) - lo) / (hi - lo), 0.0, 1.0)
    scaled = 1.0 + t * (length(IMPROVEMENT_COLORS) - 1)
    lower_idx = clamp(floor(Int, scaled), 1, length(IMPROVEMENT_COLORS))
    upper_idx = clamp(lower_idx + 1, 1, length(IMPROVEMENT_COLORS))
    frac = scaled - lower_idx
    lo_color = IMPROVEMENT_COLORS[lower_idx]
    hi_color = IMPROVEMENT_COLORS[upper_idx]
    return RGB{Float64}(
        (1.0 - frac) * red(lo_color) + frac * red(hi_color),
        (1.0 - frac) * green(lo_color) + frac * green(hi_color),
        (1.0 - frac) * blue(lo_color) + frac * blue(hi_color),
    )
end

function fig9_legend_marker_points()
    return Point2f[(0.18, 0.5), (0.50, 0.5), (0.82, 0.5)]
end

function fig9_grouped_marker_element(marker, colors; strokecolor=:transparent, strokewidth=0)
    return MarkerElement(
        marker=marker,
        color=colors,
        markersize=FIG9_LEGEND_MARKER_SIZE,
        strokecolor=strokecolor,
        strokewidth=strokewidth,
        points=fig9_legend_marker_points(),
    )
end

function add_fig9_legend!(target_slot, color_lims)
    improvement_colors = [ratio_color(value, color_lims) for value in LEGEND_IMPROVEMENT_VALUES]
    elements = [
        fig9_grouped_marker_element(:rect, reverse(CONTOUR_FILL_COLORS); strokecolor=CONTOUR_LINE_COLOR, strokewidth=1.2),
        fig9_grouped_marker_element(:circle, improvement_colors),
        fig9_grouped_marker_element(:circle, improvement_colors; strokecolor=:black, strokewidth=FIG9_MARKER_STROKE_WIDTH),
    ]
    labels = [
        "BBH Population",
        "Observed Events",
        "Inspiral SNR too low",
    ]

    return Legend(
        target_slot,
        elements,
        labels;
        framevisible=false,
        tellheight=true,
        tellwidth=false,
        halign=:left,
        valign=:top,
        patchsize=(64, 24),
        rowgap=4,
        labelsize=FIG9_LEGEND_LABELSIZE,
    )
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

function draw_top_colored_histogram!(ax::Axis, edges::Vector{Float64}, binned_ratios::Vector{Vector{Float64}}, color_lims)
    for (bin_idx, ratios) in enumerate(binned_ratios)
        for (stack_idx, ratio) in enumerate(ratios)
            y0 = Float64(stack_idx - 1)
            y1 = Float64(stack_idx)
            points = Point2f[
                (edges[bin_idx], y0),
                (edges[bin_idx + 1], y0),
                (edges[bin_idx + 1], y1),
                (edges[bin_idx], y1),
            ]
            poly!(ax, points; color=ratio_color(ratio, color_lims), strokecolor=:transparent)
        end
    end
    return ax
end

function draw_side_colored_histogram!(ax::Axis, edges::Vector{Float64}, binned_ratios::Vector{Vector{Float64}}, color_lims)
    for (bin_idx, ratios) in enumerate(binned_ratios)
        for (stack_idx, ratio) in enumerate(ratios)
            x0 = Float64(stack_idx - 1)
            x1 = Float64(stack_idx)
            points = Point2f[
                (x0, edges[bin_idx]),
                (x1, edges[bin_idx]),
                (x1, edges[bin_idx + 1]),
                (x0, edges[bin_idx + 1]),
            ]
            poly!(ax, points; color=ratio_color(ratio, color_lims), strokecolor=:transparent)
        end
    end
    return ax
end

function draw_top_histogram!(ax::Axis, base_values::Vector{Float64}, fisher_values::Vector{Float64},
    observable_values::Vector{Float64}, fisher_ratios::Vector{Float64}, observable_ratios::Vector{Float64},
    xlim::Tuple{<:Real, <:Real}, color_lims; nbins::Int=20)
    edges, base_density = histogram_bin_densities(base_values, xlim; nbins=nbins)
    _, fisher_counts = histogram_bin_counts(fisher_values, xlim; nbins=nbins)
    _, observable_counts = histogram_bin_counts(observable_values, xlim; nbins=nbins)
    fisher_binned_ratios = histogram_bin_ratios(fisher_values, fisher_ratios, edges, xlim)
    observable_binned_ratios = histogram_bin_ratios(observable_values, observable_ratios, edges, xlim)
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

    base_x, base_y = step_xy(edges, base_counts)
    lines!(ax, base_x, base_y; color=HIST_BASE_EDGE_COLOR, linewidth=1.4)
    draw_top_colored_histogram!(ax, edges, fisher_binned_ratios, color_lims)
    draw_top_colored_histogram!(ax, edges, observable_binned_ratios, color_lims)
    fisher_x, fisher_y = step_xy(edges, fisher_counts)
    observable_x, observable_y = step_xy(edges, observable_counts)
    lines!(ax, fisher_x, fisher_y; color=HIST_OUTLINE_COLOR, linewidth=HIST_LINE_WIDTH)
    lines!(ax, observable_x, observable_y; color=HIST_OUTLINE_COLOR, linewidth=HIST_OBSERVABLE_LINE_WIDTH)
    lines!(ax, [Float64(xlim[1]), Float64(xlim[2])], [0.0, 0.0]; color=:black, linewidth=1)
    return ax
end

function draw_side_histogram!(ax::Axis, base_values::Vector{Float64}, fisher_values::Vector{Float64},
    observable_values::Vector{Float64}, fisher_ratios::Vector{Float64}, observable_ratios::Vector{Float64},
    color_lims; nbins::Int=20)
    edges, base_density = histogram_bin_densities(base_values, FIG9_MC_LIMITS; nbins=nbins)
    _, fisher_counts = histogram_bin_counts(fisher_values, FIG9_MC_LIMITS; nbins=nbins)
    _, observable_counts = histogram_bin_counts(observable_values, FIG9_MC_LIMITS; nbins=nbins)
    fisher_binned_ratios = histogram_bin_ratios(fisher_values, fisher_ratios, edges, FIG9_MC_LIMITS)
    observable_binned_ratios = histogram_bin_ratios(observable_values, observable_ratios, edges, FIG9_MC_LIMITS)
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

    base_x, base_y = side_step_xy(edges, base_counts)
    lines!(ax, base_x, base_y; color=HIST_BASE_EDGE_COLOR, linewidth=1.4)
    draw_side_colored_histogram!(ax, edges, fisher_binned_ratios, color_lims)
    draw_side_colored_histogram!(ax, edges, observable_binned_ratios, color_lims)
    fisher_x, fisher_y = side_step_xy(edges, fisher_counts)
    observable_x, observable_y = side_step_xy(edges, observable_counts)
    lines!(ax, fisher_x, fisher_y; color=HIST_OUTLINE_COLOR, linewidth=HIST_LINE_WIDTH)
    lines!(ax, observable_x, observable_y; color=HIST_OUTLINE_COLOR, linewidth=HIST_OBSERVABLE_LINE_WIDTH)
    lines!(ax, [0.0, 0.0], [FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2]]; color=:black, linewidth=1)
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
    scatterplot = scatter!(ax, xvals, yvals; color=cvals, colormap=IMPROVEMENT_COLORMAP, colorrange=color_lims, markersize=FIG9_MARKER_SIZE)
    scatter!(
        ax,
        catalog[x_key][problematic],
        catalog["mc"][problematic];
        color=ratio[problematic],
        colormap=IMPROVEMENT_COLORMAP,
        colorrange=color_lims,
        markersize=FIG9_MARKER_SIZE,
        strokecolor=:black,
        strokewidth=FIG9_MARKER_STROKE_WIDTH,
    )

    xlims!(ax, Float64(xlim[1]), Float64(xlim[2]))
    ylims!(ax, FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2])
    ax.xlabel = xlabel
    ax.ylabel = ylabel
    ax.xlabelsize = FIG9_AXIS_LABELSIZE
    ax.ylabelsize = FIG9_AXIS_LABELSIZE
    ax.xticklabelsize = FIG9_AXIS_TICK_LABELSIZE
    ax.yticklabelsize = FIG9_AXIS_TICK_LABELSIZE
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

    return scatterplot
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
    max_abs_ratio = isempty(valid_ratio) ? 1.0 : maximum(abs.(valid_ratio))
    color_lims = (-max_abs_ratio, 0)

    fig = Figure(size=FIG9_SIZE, backgroundcolor=:white)
    main_width = (1.0 - fig9_side_width_fraction()) / 4.0

    top_axes = [Axis(fig[1, i], backgroundcolor=:transparent) for i in 1:4]
    colorbar_grid = GridLayout(fig[1, 5])
    main_axes = [Axis(fig[2, i], backgroundcolor=:white) for i in 1:4]
    side_ax = Axis(fig[2, 5], backgroundcolor=:white)

    for col in 1:4
        colsize!(fig.layout, col, Relative(main_width))
    end
    colsize!(fig.layout, 5, Relative(fig9_side_width_fraction()))
    rowsize!(fig.layout, 1, Relative(FIG9_TOP_HIST_FRACTION))
    rowsize!(fig.layout, 2, Relative(1.0 - FIG9_TOP_HIST_FRACTION))
    colgap!(fig.layout, FIG9_COL_GAP)
    rowgap!(fig.layout, FIG9_ROW_GAP)

    draw_top_histogram!(top_axes[1], catalog["invq"][kde_selected], catalog["invq"][fisher_hist_selected], catalog["invq"][observable_hist_selected],
        ratio[fisher_hist_selected], ratio[observable_hist_selected], (0.0, 1.0), color_lims)
    draw_top_histogram!(top_axes[2], catalog["iota"][kde_selected], catalog["iota"][fisher_hist_selected], catalog["iota"][observable_hist_selected],
        ratio[fisher_hist_selected], ratio[observable_hist_selected], (0.0, π), color_lims)
    draw_top_histogram!(top_axes[3], catalog["chi_eff"][kde_selected], catalog["chi_eff"][fisher_hist_selected], catalog["chi_eff"][observable_hist_selected],
        ratio[fisher_hist_selected], ratio[observable_hist_selected], (-1.0, 1.0), color_lims)
    draw_top_histogram!(top_axes[4], catalog["z"][kde_selected], catalog["z"][fisher_hist_selected], catalog["z"][observable_hist_selected],
        ratio[fisher_hist_selected], ratio[observable_hist_selected], (0.0, FIG9_Z_THRESHOLD), color_lims)
    foreach(ax -> hide_hist_axis!(ax; bottom_spine=true), top_axes)

    scatter_ref = build_main_panel!(main_axes[1], catalog, results, indices, "invq", (0.0, 1.0), ratio, color_lims;
        xlabel=labels["invq"], ylabel=labels["mc"], show_yticks=true, show_yticklabels=true)
    build_main_panel!(main_axes[2], catalog, results, indices, "iota", (0.0, π), ratio, color_lims;
        xlabel=labels["iota"], show_yticks=true, show_yticklabels=false)
    build_main_panel!(main_axes[3], catalog, results, indices, "chi_eff", (-1.0, 1.0), ratio, color_lims;
        xlabel=labels["chi_eff"], show_yticks=true, show_yticklabels=false)
    build_main_panel!(main_axes[4], catalog, results, indices, "z", (0.0, FIG9_Z_THRESHOLD), ratio, color_lims;
        xlabel=labels["z"], show_yticks=true, show_yticklabels=false)
    linkyaxes!(main_axes...)

    draw_side_histogram!(side_ax, catalog["mc"][kde_selected], catalog["mc"][fisher_hist_selected], catalog["mc"][observable_hist_selected],
        ratio[fisher_hist_selected], ratio[observable_hist_selected], color_lims)
    hide_hist_axis!(side_ax; left_spine=true, hide_x=true, hide_y=true)
    ylims!(side_ax, FIG9_MC_LIMITS[1], FIG9_MC_LIMITS[2])

    add_fig9_legend!(colorbar_grid[1, 1], color_lims)
    Label(
        colorbar_grid[3, 1],
        IMPROVEMENT_COLORBAR_TITLE;
        fontsize=FIG9_COLORBAR_LABELSIZE,
        halign=:center,
        valign=:bottom,
        tellwidth=false,
    )
    Colorbar(
        colorbar_grid[4, 1],
        scatter_ref;
        vertical=false,
        flipaxis=false,
        labelvisible=false,
        tickalign=0,
        ticklabelsize=FIG9_COLORBAR_TICK_LABELSIZE,
        size=18,
    )
    rowsize!(colorbar_grid, 1, Auto(0.50))
    rowsize!(colorbar_grid, 2, Fixed(FIG9_LEGEND_COLORBAR_PAD))
    rowsize!(colorbar_grid, 3, Auto(0.20))
    rowsize!(colorbar_grid, 4, Auto(0.30))
    rowgap!(colorbar_grid, 2)

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
