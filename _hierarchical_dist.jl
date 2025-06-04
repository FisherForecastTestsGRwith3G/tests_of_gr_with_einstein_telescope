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
