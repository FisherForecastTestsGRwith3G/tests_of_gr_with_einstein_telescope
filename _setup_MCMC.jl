
"""
Function that tries to simulate the hyperparamter posterior p(μ, σ | φ_k , Δ_k) with the 
help of an MCMC calculation.

DESCRIPTION:
    The function uses the Turing.jl packaged to calculate the posterior via a MCMC simulation.
    The distribution is given by the function
    
        p(μ, σ | φ_k , Δ_k) ~ p_prior(μ, σ) * Π_k exp(-(μ-dφ_k)^2/(2(Δ_k^2 + σ^2)))/sqrt(2π(Δ_k^2 + σ^2))

    The prior is defined as a flat prior i.e. a uniform distribtion with certain bounds

INPUTS:
    dphi0_k         δφ_k-values
    delta_k         Δ_k-values
    prior_mu        Prior-bounds for μ
    prior_sigma     upper prior for σ

    OPTIONAL:
    n_samples       number of samples in the chain

"""
function run_MCMC(
    dphi0_k::Vector{Float64},
    delta_k::Vector{Float64},
    prior_mu::Tuple{Float64, Float64},
    prior_sigma::Float64,
    n_samples::Int64 = 10000)

    # first check inputs
    if sum(delta_k.==0) != 0

        @warn "It seems some of the elemetns of `dphi0_k` are zero. These parameters where removed"

        delta_k = delta_k[delta_k.!=0]
        dphi0_k = dphi0_k[delta_k.!=0]
    end

    # set up bayesian model
    @model function normal_model(dphi0_k, delta_k)
        # Priors
        mu ~ Uniform(prior_mu[1], prior_mu[2])
        sigma ~ Uniform(0, prior_sigma) 
        
        # Likelihood
        for i in 1:length(dphi0_k)
            try
                dphi0_k[i] ~ Normal(mu, sqrt(sigma^2 + delta_k[i]^2))
            catch
                println("Error at ", i)
                println("mu = ", mu)
                println("sigma = ", sigma)
                println("measurement = ", dphi0_k[i])
                println("measurement error = ", delta_k[i])
            end
        end
    end

    # Perform inference
    model = normal_model(dphi0_k, delta_k)
    chain = sample(model, NUTS(), n_samples)

    return chain
end
 
"""
Plots the 2d contour of a two dimensional data set of x and y values using
a kernel density estimate.

INPUT:
    x                   x-data
    y                   y-data
    val_inj             used to mark injected values
"""
function plot2DContourKDE(x::Matrix{Float64}, y::Matrix{Float64}, val_inj::Union{Nothing, Tuple{Float64,Float64}})
    
    # Estimate the kernel density
    xy = [x y]
    kde_result = kde(xy)

    # Extract the density values
    Z = kde_result.density
    X = kde_result.x
    Y = kde_result.y

    titel_str = "MCMC"
    if val_inj != nothing
        true_mu = val_inj[1]
        true_sigma =  val_inj[2]
        titel_str = "MCMC at \$\\mu=\$$true_mu and \$\\sigma=\$$true_sigma"
    end

    MCMC_plot=plot(
        xlabel="\$\\mu\$ ",
        ylabel="\$\\sigma\$",
        title=titel_str,
        legend=:topright,  
        minorticks=9, 
        minorgrid=true, 
        grid=true, 
        minorgridwidth=2.,
        minorgridalpha=.04, 
        gridwidth=0.5, 
        gridalpha=0.5, 
        cbar=false,
        labelfontsize=21, 
        size=(800, 600),
        framestyle=:box, 
        #left_margin = 2mm, 
        #bottom_margin = 2mm, 
        #right_margin = 2.5mm, 
        #top_margin = 2.5mm,
        #TODO: This last 4 lines gave an error for me (Joachim). 
        # I commented them now because they are only cosmetic.
    )

    
    # Find the contour levels corresponding to 1σ, 2σ, and 3σ
    sorted_Z = sort(Z[:], rev=true)
    cumsum_Z = cumsum(sorted_Z) / sum(sorted_Z)  # Cumulative sum to get percentiles
    #levels = [0.682, 0.954, 0.997]  # 1σ, 2σ, 3σ intervals
    levels = [0.39, 0.86, 0.99]

    # Find density values corresponding to these probability levels
    sigma_levels = []
    for level in levels
        logical_index = findfirst(cumsum_Z .>= level)
        append!(sigma_levels, [sorted_Z[logical_index]])
    end

    # Plot contour
    contour!(MCMC_plot, X, Y, Z', levels=sigma_levels, color=:viridis, linewidth=2, label="σ levels")
    #scatter!(x, y, alpha=0.3, label="Samples")
    
    #plot injected values
    if val_inj != nothing
        true_mu = val_inj[1]
        true_sigma =  val_inj[2]
        scatter!(MCMC_plot, [true_mu], [true_sigma], label="True values")
    end
    return MCMC_plot
end

function plot2DContourKDE!(plot_, x::Matrix{Float64}, y::Matrix{Float64})

    # Estimate the kernel density
    xy = [x y]
    kde_result = kde(xy)

    # Extract the density values
    Z = kde_result.density
    X = kde_result.x
    Y = kde_result.y
    
    # Find the contour levels corresponding to 1σ, 2σ, and 3σ
    sorted_Z = sort(Z[:], rev=true)
    cumsum_Z = cumsum(sorted_Z) / sum(sorted_Z)  # Cumulative sum to get percentiles
    #levels = [0.682, 0.954, 0.997]  # 1σ, 2σ, 3σ intervals
    levels = [0.39, 0.86, 0.99]

    # Find density values corresponding to these probability levels
    sigma_levels = [sorted_Z[findfirst(cumsum_Z .>= level)] for level in levels]

    # Plot contour
    contour!(plot_, X, Y, Z', levels=sigma_levels, color=:viridis, linewidth=2, label="σ levels")
   
    # scatter!(x, y, alpha=0.3, label="Samples")
    # scatter!(MCMC_plot, [true_mu], [true_sigma], label="True values")
    return plot_
end
