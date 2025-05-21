using QuadGK
using Turing
using Distributions
using Random
using NLsolve
using SpecialFunctions

"""
Function that evaluates the posterior distribution p(δφ_n | data) for the beyond GR deformation coefficient δφ_n, given the parameters φ_k , Δ_k, for k in (1,...,N) measurements of individual GW events.

Since this is an unnormalized 1D pdf (in delta_phi), I will perform an MCMC to sample from the distribution, and then produce a violin plot with the samples. 
Care should be taken of removing the burn-in samples and make sure that the chain has converged.

"""

# Remove struct below, as it is not used
# struct Delta_phi_posterior_pdf


#     dphi0_k::Vector{Float64}
#     delta_k::Vector{Float64}

#     # Private variables

#     _a::Float64
#     _b::Float64
#     _c::Float64
#     _d::Float64

#     _normalization_prefactor::Float64
#     _sigma_maximum_value::Float64
#     _sigma_normalization_integration_range::Ntuple{2, Float64}

#     function Delta_phi_posterior_pdf(dphi0_k::Float64, delta_k::Float64)
#         new(a, b, c, d)
#     end


# end

function _evaluate_delta_phi_quantities_abcd(sigma::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64})
    # This function evaluates the quantities needed to evaluate the distribution p(δφ_n | data) for a given value of σ
    # The function returns the values of the quantities a, b, c, d which appear in the integrand

    a = sum(1.0 ./ (sigma^2 .+ delta_k.^2))
    b = sum(dphi0_k ./ (sigma^2 .+ delta_k.^2))
    c = -0.5 * sum(dphi0_k.^2 ./ (sigma^2 .+ delta_k.^2))
    d = -0.5 * sum(log1p.((sigma ./ delta_k).^2))

    return a, b, c, d
end

function _evaluate_delta_phi_sigma_integrand(sigma::Float64, delta_phi::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64})
    # This function evaluates the integrand of the distribution p(δφ_n | data), at a given point σ, given the value of δφ_n (delta_phi)

    if sigma < 0
        return 0.0
    end

    a, b, c, d = _evaluate_delta_phi_quantities_abcd(sigma, dphi0_k, delta_k)

    log_integrand = -0.5 * (a * delta_phi^2 - 2. * b * delta_phi - (b*sigma)^2)/(1 + a * sigma^2) + c + d - 0.5 * log1p(a * sigma^2)

    return exp(log_integrand)
end

function _evaluate_delta_phi_pdf(delta_phi::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64})
    # This function estimates the integral over sigma, from 0 to +infinity, of the distribution p(δφ_n | data), given the value of δφ_n (delta_phi)
    # I evaluate the integral using the Gauss-Kronrod quadrature, variables the QuadGK package

    result, err = quadgk(sigma -> _evaluate_delta_phi_sigma_integrand(sigma, delta_phi, dphi0_k, delta_k), 0., + Inf, rtol = 1e-7, atol = 0.)

    println("Integral result: ", result, " ± ", err, ", relative error: ", err / abs(result))

    return result
end

function obtain_sample_delta_phi_pdf(dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; n_samples::Int64 = 5000, burn_in::Int64 = 1000)
    # This function samples from the distribution p(δφ_n | data) for a given set of parameters φ_k , Δ_k

    struct _delta_phi_pdf{P} <: ContinuousUnivariateDistribution
        params::P
    end

    Distributions.logpdf(distr::_delta_phi_pdf, delta_phi::Real) = log(_evaluate_delta_phi_pdf(delta_phi, distr.params...))

    @model function sample_model(params)
        delta_phi ~ _delta_phi_pdf(params)
    end

    model = sample_model((dphi0_k, delta_k))

    #Use the Metropolis-Hastings algorithm to sample from the posterior distribution
    chain = sample(model, NUTS(), n_samples + burn_in)

    #Remove the burn-in samples
    samples = chain[:delta_phi][burnin+1:end]

    return samples
end


"""
Function that calculates the posterior distribustion p(delta_phi | data, σ=0) for the beyond GR deformation coefficient δφ_n, given the parameters φ_k , Δ_k, for k in (1,...,N) measurements of individual GW events, 
starting from the hierarchical hyperparameter distribution p(μ, σ | φ_k , Δ_k), but in this case, conditioned with σ=0.

DESCRIPTION:
For an event with index k we have
    φ_k ... the peak/mean of the distribution
    Δ_k ... the spread/standard deviation of the distribution

Here we assume that σ=0, i.e. that the deviation in GR parameters must be the same in every observation,

INPUT:
    dphi0_k    φ_k-values
    delta_k    Δ_k-values

"""
function getMuDistTIGER(dphi0_k::Vector{Float64}, delta_k::Vector{Float64})
    
    a0 = sum(1.0 ./ delta_k.^2)
    b0 = sum(dphi0_k ./ delta_k.^2)
    #c0 = -0.5 * sum(dphi0_k.^2 ./ delta_k.^2)

    meanMuDist = b0/a0
    stdMuDist = sqrt(1/a0)

    #The pdf for mu is given by Normal distribution with mean and std calculated above
    #log_p_mu = -0.5*((mu - meanMuDist)^2)/(stdMuDist^2) - 0.5*log(2*pi*stdMuDist^2)
    # i will return the mean and variance of this distribution in full generality
    return meanMuDist, stdMuDist
end

function getUpperLimitOnAbsValueNormalVariable(mean_dst, std_dev_dst, probability_value = 0.9)
    
    # I assume a gaussian distribution for x, with mean and std given by the input arguments
    # so p(x | mean_dst, std_dst) = 1/(sqrt(2*pi)*std_dst)*exp(-0.5*((x - mean_dst)/std_dst)^2)

    # This function will return the (positive) value of x_upper_limit, such that the probability of the absolute value of x being less than x_upper_limit is probability_value (e.g probability_value = 0.9)
    # i.e. P(|x| < x_upper_limit) = probability_value

    # In practice I want to calculate the value of x_upper_limit such that the integral from -x_upper_limit to x_upper_limit of the pdf is probability_value
    # i.e. probability_value = integrate(1/(sqrt(2*pi)*std_dst)*exp(-0.5*((x - mean_dst)/std_dst)^2), x, -xLim, xLim)

    functionToSolve(xLim) = 1/2. * (erf.((xLim .- mean_dst)/(sqrt(2) * std_dev_dst)) + erf.((xLim .+ mean_dst)/(sqrt(2)  *std_dev_dst))) .- probability_value

    initial_guess = [(abs(mean_dst) + std_dev_dst)]
    solution = nlsolve(functionToSolve, initial_guess)
    x_upper_limit = solution.zero[1]

    return x_upper_limit
end

"""
DESCRIPTION:
Function that calculates the 90% upper limit of the beyond GR deformation coefficient δφ_n, given the parameters φ_k , Δ_k, for k in (1,...,N) measurements of individual GW events. 

The function evaluates the 90% upper limit by evauating the probability p(δφ_n | data), and then evaluating the value of δφ_n_limit such that the probability of δφ_n being less than δφ_n_limit is 0.9.

The probability p(δφ_n | data) is assumed to come from the hierarchical hyperparamter distribution pHier(μ, σ | φ_k , Δ_k), conditioned with σ=0.
From p(δφ_n | data) = ∫ p(δφ_n | μ, σ) pHier(μ, σ | data) dμ dσ it follows that p(δφ_n | data) = pHier(μδφ_n, 0 | data) dμ dσ .

Therefore p(δφ_n | data) is given by a gaussian distribution, with a mean different from zero. To evaluate the quantity of interest δφ_n_limit, such that p(|δφ_n| < δφ_n_limit| data) < 0.9, then we evaluate the integral
∫_(-δφ_n_limit)^(δφ_n_limit) p(δφ_n | data) dδφ_n = 0.9.


INPUT:
    dphi0_k    φ_k-values
    delta_k    Δ_k-values

    Optional
    upper_limit_probability_value   Value of the probability that the value of δφ_n is less than δφ_n_limit.
               value 0.9 by default. 
"""
function get90PctUpperLimitMuDistTIGER(dphi0_k::Vector{Float64}, delta_k::Vector{Float64}, upper_limit_probability_value = 0.9)
    # Returns the 90% upper limit (or upper_limit_probability_value) for the "naive" TIGER approach (actually the one where we set sigma = 0, and so impose all events to have the same beyond GR deviation)
    # takes as input the dphi0_k and delta_k values for the events, at a given PN order
    
    meanMuDist, stdMuDist  = getMuDistTIGER(dphi0_k, delta_k)
    upper_limit = getUpperLimitOnAbsValueNormalVariable(meanMuDist, stdMuDist, upper_limit_probability_value)

    return upper_limit
end

