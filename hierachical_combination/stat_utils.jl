"""
Determine one or several quantiles from a one-dimensional posterior known on a
grid (OGD = OnGridDist).

The values in `p_values` are interpreted as probability masses on the grid
`x_values`. If `weights` is provided, it is multiplied into `p_values` before
constructing the cumulative distribution. This is useful when the grid is
non-uniform and the bin widths should be accounted for explicitly.

The method with `quantile_levels::AbstractVector` computes the cumulative
distribution once and returns a vector of quantiles in the same order as the
requested levels. The scalar method is a thin wrapper around the vector
implementation.
"""
function quantile1dOGD(
    x_values::AbstractVector{<:Real},
    p_values::AbstractVector{<:Real},
    quantile_levels::AbstractVector{<:Real};
    weights=nothing,
    )

    length(x_values) == length(p_values) ||
        throw(ArgumentError("`x_values` and `p_values` must have the same length."))
    length(x_values) >= 2 ||
        throw(ArgumentError("At least two grid points are required."))

    x_grid = Float64.(x_values)
    p_grid = Float64.(p_values)
    q_grid = Float64.(quantile_levels)

    issorted(x_grid) || throw(ArgumentError("`x_values` must be sorted in increasing order."))
    any(p_grid .< 0.0) && throw(ArgumentError("`p_values` must be non-negative."))
    (any(q_grid .< 0.0) || any(q_grid .> 1.0)) &&
        throw(ArgumentError("All `quantile_levels` must lie between 0 and 1."))

    if !isnothing(weights)
        length(weights) == length(x_grid) ||
            throw(ArgumentError("`weights` must have the same length as `x_values`."))
        weights = Float64.(weights)
    end
    !isnothing(weights) && any(weights .< 0.0) &&
        throw(ArgumentError("`weights` must be non-negative."))

    probability_masses = if isnothing(weights)
        p_grid
    else
        p_grid .* weights
    end
    total_probability = sum(probability_masses)
    total_probability > 0.0 || throw(ArgumentError("`p_values` must contain positive probability mass."))

    cdf = cumsum(probability_masses) ./ total_probability
    cdf ./= cdf[end]

    quantiles = zeros(length(q_grid))
    for (idx_q, quantile_level) in enumerate(q_grid)
        if quantile_level == 0.0
            quantiles[idx_q] = x_grid[1]
            continue
        elseif quantile_level == 1.0
            quantiles[idx_q] = x_grid[end]
            continue
        end

        idx = searchsortedfirst(cdf, quantile_level)
        idx = clamp(idx, 2, length(cdf))

        x_left = x_grid[idx - 1]
        x_right = x_grid[idx]
        cdf_left = cdf[idx - 1]
        cdf_right = cdf[idx]

        if cdf_right == cdf_left
            quantiles[idx_q] = x_right
            continue
        end

        frac = (quantile_level - cdf_left) / (cdf_right - cdf_left)
        quantiles[idx_q] = x_left + frac * (x_right - x_left)
    end

    return quantiles
end

function quantile1dOGD(
    x_values::AbstractVector{<:Real},
    p_values::AbstractVector{<:Real},
    quantile_level::Float64;
    weights=nothing,
    )

    return only(quantile1dOGD(x_values, p_values, [quantile_level]; weights=weights))
end

"""
Calculate the posterior level corresponding to a requested credible interval on a
discrete two-dimensional distribution given on a grid (OGD = OnGridDist).

The grid values in `p_values` are treated as bin densities. The optional
`weights` argument can be used to supply bin areas for non-uniform grids.
"""
function getContourLevelOGD(
    p_values::AbstractArray{<:Real},
    CI::Float64=0.90;
    weights=nothing,
    )

    0.0 < CI < 1.0 || throw(ArgumentError("`CI` must be strictly between 0 and 1."))

    density_values = vec(Float64.(p_values))
    any(density_values .< 0.0) &&
        throw(ArgumentError("`p_values` must be non-negative."))

    probability_masses = if isnothing(weights)
        density_values
    else
        size(weights) == size(p_values) ||
            throw(ArgumentError("`weights` must have the same size as `p_values`."))
        weight_values = vec(Float64.(weights))
        any(weight_values .< 0.0) &&
            throw(ArgumentError("`weights` must be non-negative."))
        density_values .* weight_values
    end

    total_probability = sum(probability_masses)
    total_probability > 0.0 || throw(ArgumentError("`p_values` must contain positive probability mass."))

    sorted_indices = sortperm(density_values; rev=true)
    cumulative_probability = cumsum(probability_masses[sorted_indices]) ./ total_probability
    level_idx = findfirst(cumulative_probability .>= CI)

    isnothing(level_idx) &&
        throw(ArgumentError("Could not determine the contour level for the requested `CI`."))

    return density_values[sorted_indices[level_idx]]
end

"""
Calculate the enclosed probability of the iso-volume region defined by
`p_values >= p_ref` for a distribution known on a grid (OGD = OnGridDist).

The values in `p_values` are interpreted as grid-based probability densities or
probability masses. If `weights` is provided, it is multiplied into `p_values`
before summing, so the enclosed probability can account for non-uniform
grid-cell volumes.

The returned value is the fraction of the total probability mass contained in
all cells whose value is at least `p_ref`.
"""
function getEnclosedIsoVolProbOGD(
    p_values::AbstractArray{<:Real},
    p_ref::Float64,
    ; weights=nothing,
)

    bin_values = if isnothing(weights)
        vec(Float64.(p_values))
    else
        size(weights) == size(p_values) ||
            throw(ArgumentError("`weights` must have the same size as `p_values`."))
        vec(Float64.(p_values) .* Float64.(weights))
    end

    total_probability = sum(bin_values)
    total_probability > 0.0 || throw(ArgumentError("`p_values` must contain positive probability mass."))

    enclosed_bins = vec(Float64.(p_values) .>= p_ref)
    enclose_probability_mass = sum(bin_values[enclosed_bins])

    enclosed_probability = enclose_probability_mass ./ total_probability

    return enclosed_probability
end
