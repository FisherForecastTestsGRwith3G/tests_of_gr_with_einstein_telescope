"""
Determine effective bin widths for a one-dimensional, strictly increasing grid.
"""
function get1dGridWeights(x_values::AbstractVector{<:Real})

    length(x_values) >= 2 ||
        throw(ArgumentError("At least two grid points are required."))

    x_grid = Float64.(x_values)
    issorted(x_grid) || throw(ArgumentError("`x_values` must be sorted in increasing order."))

    widths = zeros(length(x_grid))
    widths[1] = 0.5 * (x_grid[2] - x_grid[1])
    widths[end] = 0.5 * (x_grid[end] - x_grid[end - 1])

    for idx in 2:length(x_grid) - 1
        widths[idx] = 0.5 * (x_grid[idx + 1] - x_grid[idx - 1])
    end

    any(widths .<= 0.0) && throw(ArgumentError("`x_values` must be strictly increasing."))

    return widths
end

"""
Determine effective cell areas for a two-dimensional, strictly increasing grid.
"""
function get2dGridWeights(
    x_values::AbstractVector{<:Real},
    y_values::AbstractVector{<:Real},
    )

    x_weights = get1dGridWeights(x_values)
    y_weights = get1dGridWeights(y_values)

    return x_weights .* transpose(y_weights)
end

"""
    buildCoarseHyperparamDomain(hyper, center_mu=0.0, center_sig=0.0, ks_hyper=5)

Construct a coarse hyperparameter grid for the hierarchical mean `mu` and width
`sigma` from a `hyperparamDistTIGER` object.

The grid spacing is derived from `hyper.delta_k`. The returned `mu` grid is
approximately symmetric around `center_mu`, while the returned `sigma` grid is
shifted by `center_sig` and truncated so that it only contains non-negative
values.

Returns `(mu, sigma)`.
"""
function buildCoarseHyperparamDomain(
    hyper::hyperparamDistTIGER,
    center_mu::Float64 = 0.0,
    center_sig::Float64 = 0.0,
    ks_hyper::Int = 5, 
    )

    delta_k = hyper.delta_k

    spread     = maximum(delta_k)  / sqrt(length(delta_k))
    mu_min     = center_mu  - ks_hyper * spread
    mu_max     = center_mu  + ks_hyper * spread
    # IMPROVEME: This is not a very ideal estimate for sigma spread
    # It is solely based on the idea that the distribution fits kind 
    # of in a square. Typically the sigma_direction is spread out more.  
    sigma_min  = center_sig - ks_hyper * spread
    sigma_max  = center_sig + ks_hyper * spread

    shell_point_counts = (200, 100, 50, 25, 25)

    function build_positive_shell_grid(max_shell::Int)
        grid = [0.0]
        for idx_shell in 1:max_shell
            shell_start = (idx_shell - 1) * spread
            shell_stop = idx_shell * spread
            n_shell_points = shell_point_counts[min(idx_shell,5)]
            shell_grid = collect(range(shell_start, shell_stop; length=n_shell_points + 1))
            append!(grid, shell_grid[2:end])
        end
        return grid 
    end

    mu_pos = build_positive_shell_grid(ks_hyper)
    mu_neg = -reverse(mu_pos[2:end])
    mu = vcat(mu_neg, mu_pos) .+ center_mu
    sigma_pos = build_positive_shell_grid(ks_hyper)
    sigma_neg = Float64[]
    sigma = vcat(sigma_neg, sigma_pos) .+ center_sig
    sigma = sigma[searchsortedfirst(sigma, 0.0):end]

    return mu, sigma
end

"""
    buildUniform1dGridWithCenter(x_min, x_max, n_points, x0=0.0)

Build a one-dimensional uniform grid with `n_points` between `x_min` and
`x_max`, ensuring that `0.0` is included exactly. The resulting grid is then
shifted by `x0`.

Returns a vector of grid points.
"""
function buildUniform1dGridWithCenter(x_min::Float64, x_max::Float64, n_points::Int, x0::Float64=0.0)
    x_min <= 0.0 <= x_max || throw(ArgumentError("The refined grid must contain 0.0."))
    n_points >= 3 || throw(ArgumentError("The refined grid needs at least 3 points."))

    if x_min == 0.0
        return collect(range(0.0, x_max; length=n_points))
    elseif x_max == 0.0
        return collect(range(x_min, 0.0; length=n_points))
    end

    n_left = clamp(round(Int, (n_points - 1) * (-x_min) / (x_max - x_min)), 1, n_points - 2)
    n_right = n_points - n_left - 1

    left_grid = collect(range(x_min, 0.0; length=n_left + 1))
    right_grid = collect(range(0.0, x_max; length=n_right + 1)) 

    return vcat(left_grid[1:end-1], right_grid) .+ x0
end

"""
Check whether a `(mu, sigma)` posterior grid is sufficiently extended.

The diagnostic compares the full normalization of `p_mu_sigma` against the
normalization obtained after trimming away:
- 2 layers on the left in the `mu` direction
- 2 layers on the right in the `mu` direction
- 2 layers at the top in the `sigma` direction

If the relative change in normalization is smaller than `threshold_percent`,
the grid is considered good.
"""
function checkPosteriorGridSuitability(
    mu::AbstractVector{<:Real},
    sigma::AbstractVector{<:Real},
    p_mu_sigma::AbstractMatrix{<:Real};
    threshold_percent::Float64=0.001,
    )

    size(p_mu_sigma) == (length(mu), length(sigma)) ||
        throw(ArgumentError("`p_mu_sigma` must have size `(length(mu), length(sigma))`."))
    length(mu) >= 5 ||
        throw(ArgumentError("The mu-grid must contain at least 5 points to trim two layers on each side."))
    length(sigma) >= 3 ||
        throw(ArgumentError("The sigma-grid must contain at least 3 points to trim two layers at the top."))
    threshold_percent >= 0.0 ||
        throw(ArgumentError("`threshold_percent` must be non-negative."))

    full_norm = trapz((mu, sigma), p_mu_sigma)
    full_norm > 0.0 || throw(ArgumentError("`p_mu_sigma` must contain positive probability mass."))

    mu_trim = mu[3:end-2]
    sigma_trim = sigma[1:end-2]
    p_trim = p_mu_sigma[3:end-2, 1:end-2]
    trimmed_norm = trapz((mu_trim, sigma_trim), p_trim)

    norm_change = abs(full_norm - trimmed_norm)
    norm_change_percent = 100.0 * norm_change / full_norm
    grid_is_good = norm_change_percent < threshold_percent

    return (
        grid_is_good=grid_is_good,
        threshold_percent=threshold_percent,
        full_norm=full_norm,
        trimmed_norm=trimmed_norm,
        norm_change=norm_change,
        norm_change_percent=norm_change_percent,
    )
end

"""
    findOptimalGrid(hyper, center_mu=0.0, center_sig=0.0)

Construct a coarse `(mu, sigma)` grid for `hyper`, enlarge it until the
posterior support is sufficiently contained, and extract a bounding box around
the 99.99% credible region.

Returns `(mu_box_min, mu_box_max, sigma_box_min, sigma_box_max)`.
"""
function findOptimalGrid(
    hyper::hyperparamDistTIGER,
    center_mu::Float64 = 0.0,
    center_sig::Float64 = 0.0;
    verbose::Bool=true
    )

    ks_hyper = 5
    grid_check = nothing
    mu_grid = Float64[]
    sigma_grid = Float64[]
    p_mu_sigma = Matrix{Float64}(undef, 0, 0)

    while true
        mu_grid, sigma_grid = buildCoarseHyperparamDomain(
            hyper, 
            center_mu,
            center_sig,
            ks_hyper
        )

        p_mu_sigma, _, _, _, _ = getDistributionOnGrid(
            mu_grid,
            sigma_grid,
            hyper,
        )

        grid_check = checkPosteriorGridSuitability(mu_grid, sigma_grid, p_mu_sigma)
        grid_check.grid_is_good && break

        verbose && println("Increasing ks_hyper to $(ks_hyper + 1).")
        ks_hyper += 1
        ks_hyper > 10 && throw(ErrorException("Coarse hyperparameter grid search did not converge by ks_hyper = 10."))
    end

    coarse_grid_weights = get2dGridWeights(mu_grid, sigma_grid)
    ci_9999_level = getContourLevelOGD(p_mu_sigma, 0.9999; weights=coarse_grid_weights)
    ci_9999_region = findall(>=(ci_9999_level), p_mu_sigma)
    isempty(ci_9999_region) && throw(ErrorException("Failed to locate the 99.99% credible region on the coarse grid."))

    mu_idx_min = first(ci_9999_region)[1]
    mu_idx_max = mu_idx_min
    sigma_idx_min = first(ci_9999_region)[2]
    sigma_idx_max = sigma_idx_min

    for idx in ci_9999_region
        mu_idx = idx[1]
        sigma_idx = idx[2]
        mu_idx < mu_idx_min && (mu_idx_min = mu_idx)
        mu_idx > mu_idx_max && (mu_idx_max = mu_idx)
        sigma_idx < sigma_idx_min && (sigma_idx_min = sigma_idx)
        sigma_idx > sigma_idx_max && (sigma_idx_max = sigma_idx)
    end

    mu_box_min = min(mu_grid[mu_idx_min], 0.0)
    mu_box_max = max(mu_grid[mu_idx_max], 0.0)
    sigma_box_min = min(sigma_grid[sigma_idx_min], 0.0)
    sigma_box_max = max(sigma_grid[sigma_idx_max], 0.0)

    if mu_box_min == mu_box_max
        mu_box_min = mu_grid[max(mu_idx_min - 1, 1)]
        mu_box_max = mu_grid[min(mu_idx_max + 1, length(mu_grid))]
    end
    if sigma_box_min == sigma_box_max
        sigma_box_min = sigma_grid[max(sigma_idx_min - 1, 1)]
        sigma_box_max = sigma_grid[min(sigma_idx_max + 1, length(sigma_grid))]
    end
    if sigma_box_max == 0.0
        sigma_box_max = sigma_grid[min(sigma_idx_max + 1, length(sigma_grid))]
    end

    return mu_box_min, mu_box_max, sigma_box_min, sigma_box_max
end
