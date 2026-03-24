
"""
Given a Gaussian distribution N(x, mu, sigma), this function determines
the value x90 > 0, such that

    int_{-x90}^{x90} dx N(x, mu, sigma) = CI = 0.9
"""
function estimate0Symmetric90CiGaussian(mu::Float64, sigma::Float64, CI::Float64 = 0.9)

    p = Normal(mu, sigma)
    f(x) = cdf(p, x) - cdf(p, -x) - 0.9
    x90 = find_zero(f, (0.0, abs(mu) + 10sigma))
    return x90
end 

function estimate0Symmetric90CiGaussian(
    mu::Vector{Float64},
    sigma::Vector{Float64},
    CI::Float64 = 0.9,
    )

    length(mu) == length(sigma) || throw(ArgumentError("`mu` and `sigma` must have the same length."))
    return estimate0Symmetric90CiGaussian.(mu, sigma, CI)
end

"""
Calculate the effective mean value and effective standard deviation, describing 
a product of Gaussian distributions as a Gaussian distribution i.e.

    Prod_i N(x,mu_i,sig_i) = C*N(x, mu_eff, sig_eff) ,

where C is a normalization constant.
"""
function muStdEff4ProdNormal(mu::Vector{Float64}, sig::Vector{Float64})

    sig_m2 = 1.0 ./ sig.^2
    a = sum(sig_m2)
    b = sum(mu .* sig_m2)
   
    mu_eff  = b / a
    sig_eff = sqrt(1.0 / a)
    
    return mu_eff, sig_eff
end 
