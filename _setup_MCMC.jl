#using LinearAlgebra

function run_MCMC(parameters, errors, true_mu, true_sigma, chain_points)

    println("mu = ", true_mu, " sigma = ", true_sigma)


    measurement_errors = errors[errors.!=0]
    par_non_zeros = parameters[errors.!=0]
    # measurement_errors_shifted = map(M -> M * sqrt(inv(M)), measurement_errors .* measurement_errors)

    N = length(par_non_zeros)


    # measurements = par_non_zeros .+ measurement_errors_shifted .* randn(N)
    measurements = par_non_zeros .+ measurement_errors .* randn(N)

    # Bayesian model
    @model function normal_model(measurements, measurement_errors)
        # Priors
        mu ~ Normal(0, 10*true_sigma)   
        sigma ~ truncated(Normal(0, 5), 0, Inf)  # Positive sigma
        
        # Likelihood
        for i in 1:length(measurements)
            try
                measurements[i] ~ Normal(mu, sqrt(sigma^2 + measurement_errors[i]^2))
            catch
                println("Error at ", i)
                println("mu = ", mu)
                println("sigma = ", sigma)
                println("measurement = ", measurements[i])
                println("measurement error = ", measurement_errors[i])
            end
        end
    end



    # Perform inference
    model = normal_model(measurements, measurement_errors)

    #if MCMC_multiple_parallel_chains
    #    println("Running ", MCMC_n_chains_multiple_chains, " multiple chains in parallel")
    #    chain = sample(model, NUTS(), MCMCThreads(), chain_points, MCMC_n_chains_multiple_chains; init_params=[true_mu, true_sigma+1e-7]) # + 1e-7 to avoid sigma = 0
    #else
        chain = sample(model, NUTS(), chain_points, init_params = [true_mu, true_sigma+1e-7]) # + 1e-7 to avoid sigma = 0
    #end

    return chain
end
  

function plot_2d_contour(x, y, true_mu, true_sigma)

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
    #println(cumsum_Z)

    #levels = [0.682, 0.954, 0.997]  # 1σ, 2σ, 3σ intervals
    levels = [0.39, 0.86, 0.99]

    # Find density values corresponding to these probability levels
    sigma_levels = []
    for level in levels
        logical_index = findfirst(cumsum_Z .>= level)
        print("index = $(logical_index)")
        append!(sigma_levels, [sorted_Z[logical_index]])
    end
    println(sigma_levels)

    # Plot contour
    cc=contour(X, Y, Z', levels=sigma_levels, color=:viridis, linewidth=2, label="σ levels")
    #scatter!(x, y, alpha=0.3, label="Samples")
    xlabel!(L"\mu")
    ylabel!(L"\sigma")
    title!("1σ, 2σ, 3σ Intervals posteriors for μ = $(round(true_mu, sigdigits=4)) and σ = $(round(true_sigma, sigdigits=4))")
    scatter!([true_mu], [true_sigma], label="True values")
    plot!(size=(800, 600))
    return cc
end
