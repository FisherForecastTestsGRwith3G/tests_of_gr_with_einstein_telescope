function run_MCMC(parameters, errors, true_mu, true_sigma)

    println("mu = ", true_mu, " sigma = ", true_sigma)


    measurement_errors = errors[errors.!=0]
    par_non_zeros = parameters[errors.!=0]

    N = length(par_non_zeros)

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


    chain = sample(model, NUTS(), 2000, init_params = [true_mu, true_sigma+1e-7]) # + 1e-7 to avoid sigma = 0
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

    MCMC_plot=plot(xlabel="\$\\mu\$ ", ylabel="\$\\sigma\$", title="MCMC at \$\\mu=\$$true_mu and \$\\sigma=\$$true_sigma",
    legend=:topright,  minorticks=9, minorgrid=true, grid=true, minorgridwidth=2.,
    minorgridalpha=.04, gridwidth=0.5, gridalpha=0.5, cbar=false,
    labelfontsize=21, size=(800, 600), left_margin = 2mm, bottom_margin = 2mm, right_margin = 2.5mm, top_margin = 2.5mm,
    framestyle=:box)

    
    # Find the contour levels corresponding to 1σ, 2σ, and 3σ
    sorted_Z = sort(Z[:], rev=true)
    cumsum_Z = cumsum(sorted_Z) / sum(sorted_Z)  # Cumulative sum to get percentiles
    println(cumsum_Z)
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
    contour!(MCMC_plot, X, Y, Z', levels=sigma_levels, color=:viridis, linewidth=2, label="σ levels")
    #scatter!(x, y, alpha=0.3, label="Samples")
    scatter!(MCMC_plot, [true_mu], [true_sigma], label="True values")
    return MCMC_plot
end




function plot_2d_contour(x, y, true_mu, true_sigma, plot_)

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
    #scatter!(x, y, alpha=0.3, label="Samples")
    # scatter!(MCMC_plot, [true_mu], [true_sigma], label="True values")
    return plot_
end
