"""
Extract contour-line segments from a hyperparameter distribution evaluated on a
`(mu, sigma)` grid.

Returns two vectors of vectors containing the `mu` and `sigma` coordinates of
each contour segment at the requested credible level.
"""
function getCredibleContourOGD(
    mu_values::AbstractVector{<:Real},
    sigma_values::AbstractVector{<:Real},
    p_mu_sigma::AbstractMatrix{<:Real},
    CI::Float64=0.90;
    weights=nothing,
    )

    size(p_mu_sigma) == (length(mu_values), length(sigma_values)) ||
        throw(ArgumentError("`p_mu_sigma` must have size `(length(mu_values), length(sigma_values))`."))

    contour_level = getContourLevelOGD(p_mu_sigma, CI; weights=weights)

    contour_lines_mu = Vector{Vector{Float64}}()
    contour_lines_sigma = Vector{Vector{Float64}}()

    for contour_level_set in con.levels(con.contours(mu_values, sigma_values, p_mu_sigma, [contour_level]))
        for line in con.lines(contour_level_set)
            mu_line, sigma_line = con.coordinates(line)
            push!(contour_lines_mu, Float64.(mu_line))
            push!(contour_lines_sigma, Float64.(sigma_line))
        end
    end

    return contour_lines_mu, contour_lines_sigma
end
