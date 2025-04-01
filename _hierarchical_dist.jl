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
Obsolete. Please delete
"""
function OldHyperparamDistTIGER(mu::Vector{Float64}, sigma::Vector{Float64}, dphi_k::Vector{Float64}, delta_k::Vector{Float64})

    # check inputs 
    n_mu = length(mu)
    n_sigma = length(sigma)
    
    n_data = length(dphi_k)
    if length(delta_k) != n_data
        print("Problem")
    end

    # reshape data
    dphi_k = dphi_k[:]
    delta_k = delta_k[:]

    # initiate distribution arrays
    p_mu_sig = Array{Float64}(undef, (n_mu, n_sigma))
    p_sig = Array{Float64}(undef, n_sigma)

    # calculate the scale of the distribution, marginalized over mu and at z = 0
    # this scales with the number of samples i.e. n_data and cancles the growing of 
    # the exponent later such that things remain numerically stable for abitrary 
    # number n_data
    a0 = sum(1.0 ./ delta_k.^2)
    unity = ones(size(delta_k))

    b0 = sum(dphi_k ./ delta_k.^2)
    c0 = -0.5 * sum(dphi_k.^2 ./ delta_k.^2)
    f0 = 0.5 * b0^2 / a0 + c0  
    g0 = a0 / (2.0 .* pi) 
    ln_p0 =  f0 - log(g0)/2

    # calculate the actual distribution over the 2d plane
    for idx_s in 1:n_sigma
    
        sig = sigma[idx_s]
        
        if sig > 0.0 
                        
            ak = 1.0 ./ (1.0 ./ delta_k.^2 .+ 1.0 ./ sig.^2)
            aks = ak ./ sig^2 
            a = sum( 1.0 .- ak ./ sig.^2) / sig^2 
            b = sum( ak ./ delta_k.^2 .* dphi_k) / sig^2 
            c = 0.5 * sum( (ak./ delta_k.^2 .- 1.0) .* dphi_k.^2 ./ delta_k.^2) 
            
            f = -0.5 .* a .* mu.^2 .+ b .* mu .+ c 
            ln_g = 0.5 * sum(log.(aks)) 
            
            p_mu_sig[:, idx_s] = exp.(f .+ ln_g .- ln_p0)
            p_sig[idx_s] = exp(0.5*b^2/a + c .+ ln_g .- ln_p0 - 0.5 .* log.(a/(2*pi)))

        # use a different form of the distribution, since above form is ill-defined if sig == 0
        elseif  sig == 0.0 

            f = -0.5 .* a0 .* mu.^2 .+ b0 .* mu .+ c0
            p_mu_sig[:, idx_s] = exp.(f .- ln_p0)
            p_sig[idx_s] = 1

        else 
            throw(ArgumentError("Argument 'sigma' must contain only positive values!"))
        end
    end 

    # marginalize sigma by numerical integration
    p_mu = Array{Float64}(undef, n_mu)
    for idx in 1:n_mu
        p_mu[idx] = trapz(sigma, p_mu_sig[idx, :])
    end
    
    # now we need to calculate the actual normalization by integrating out mu
    n = trapz(sigma, p_sig)
    nn = trapz((mu, sigma), p_mu_sig)

    return p_mu_sig, p_sig, p_mu, n, nn 
    
end
