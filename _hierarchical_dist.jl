using NLsolve
using SpecialFunctions

"""
Function that calculates the values of the LSA hyperparameter distribution p(μ, σ | φ_k , Δ_k)
given the k in (1,...,N) measurements of individual GW events. 

DESCRIPTION:
The single event measurements in the LSA are gaussians and thus characterized by only two parameter
For an event with index k we have
    φ_k ... the peak/mean of the distribution
    Δ_k ... the spread/standard deviation of the distribution

INPUT:
    mu         μ-values
    sigma      σ-values 
    dphi0_k    φ_k-values
    delta_k    Δ_k-values
"""
function hyperparamDistTIGER(mu::Vector{Float64}, sigma::Vector{Float64}, dphi0_k::Vector{Float64}, delta_k::Vector{Float64})

    n_mu = length(mu)
    n_sigma = length(sigma)
    
    #initialize pdfs 
    log_p_mu_sigma = zeros(n_mu, n_sigma)  # 2d distribution, dim=1 ~ mu, dim=2 ~ sigma 
    log_p_sigma = zeros(n_sigma)           # 1d marginal distribution for sigma

    a0 = sum(1.0 ./ delta_k.^2)
    b0 = sum(dphi0_k ./ delta_k.^2)
    c0 = -0.5 * sum(dphi0_k.^2 ./ delta_k.^2)
    ln_a0 = log(a0) - log(2*pi)
    ln_norm_from_0 = 0.5*b0^2/a0 + c0 - 0.5*ln_a0 

    # calculate the distribution
    for idx_sigma in 1:n_sigma
        
        sig = sigma[idx_sigma]
        
        if sig >= 0

            denom = sig.^2 .+ delta_k.^2 
            a = sum(1.0 ./ denom) 
            b = sum(dphi0_k ./ denom)
            c = -0.5 * sum(dphi0_k.^2 ./ denom)
            ln_norm = 0.5*sum(log.(delta_k.^2 ./ denom))
            ln_a = log(a) - log(2*pi)

            log_p_mu_sigma[:, idx_sigma] = -0.5*mu.^2*a .+ b*mu .+ c .+ ln_norm .- ln_norm_from_0
            log_p_sigma[idx_sigma] = 0.5*b^2/a .+ c .- 0.5*ln_a .+ ln_norm .- ln_norm_from_0
        
        else

            throw(ArgumentError("Argument 'sigma' must contain only positive values!"))
        end

    end

    p_mu_sigma =  exp.(log_p_mu_sigma)
    p_sigma = exp.(log_p_sigma)

    # marginalize sigma by numerical integration
    p_mu = zeros(n_mu)
    for idx_mu in 1:n_mu
        p_mu[idx_mu] = trapz(sigma, p_mu_sigma[idx_mu,:])
    end
     
    # now we need to calculate the actual normalization by integrating out mu
    n = trapz(sigma, p_sigma)
    nn = trapz((mu, sigma), p_mu_sigma)  
                            
    return p_mu_sigma, p_sigma, p_mu, n, nn
end

"""
Function that calculates the values of the LSA hyperparameter distribution p(μ | σ=0, φ_k , Δ_k)
given the k in (1,...,N) measurements of individual GW events. 

DESCRIPTION:
The single event measurements in the LSA are gaussians and thus characterized by only two parameter
For an event with index k we have
    φ_k ... the peak/mean of the distribution
    Δ_k ... the spread/standard deviation of the distribution

Assuming σ=0, i.e. that the deviation in GR parameters must be the same in every observation, one 
obtains a different distribution, then when relaxing the assumption of having σ an unknown. This distribution is contained in the function 'hyperparamDistTIGER' at the value σ=0. However, we provide here a separate implementation.

INPUT:
    mu         μ-values
    dphi0_k    φ_k-values
    delta_k    Δ_k-values

    Optional
    norm_max   Value of the maximum of the distribution. The distribution is normalized such that the maximum 
               value 1 by default. 
"""
function naiveMuDistTIGER(mu::Vector{Float64}, dphi0_k::Vector{Float64}, delta_k::Vector{Float64}, norm_max::Float64=1.)
    
    a0 = sum(1.0 ./ delta_k.^2)
    b0 = sum(dphi0_k ./ delta_k.^2)
    c0 = -0.5 * sum(dphi0_k.^2 ./ delta_k.^2)

    log_p_mu =  -0.5*mu.^2*a0 .+ b0*mu .+ c0
    max_log_p_mu = maximum(log_p_mu)

    return exp.(log_p_mu .- max_log_p_mu)*norm_max
end

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
