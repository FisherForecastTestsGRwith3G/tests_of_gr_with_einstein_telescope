using QuadGK
using Turing
using Distributions
using Random
using NLsolve
using SpecialFunctions

# # Definitions needed for Turing MCMC in function obtain_samples_delta_phi_pdf
# struct _delta_phi_pdf{P} <: ContinuousUnivariateDistribution
#     params::P
# end

# Distributions.logpdf(distr::_delta_phi_pdf, delta_phi::Real) = log(_evaluate_delta_phi_pdf(delta_phi, distr.params...))

# Lock for plotting, needed to avoid multiple threads trying to plot at the same time (plot backend is not thread-safe)
const plot_lock = ReentrantLock()

function _evaluate_delta_phi_quantities_abcd(sigma::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64})
    # This function evaluates the quantities needed to evaluate the distribution p(δφ_n | data) for a given value of σ
    # The function returns the values of the quantities a, b, c, d which appear in the integrand

    a = sum(1.0 ./ (sigma^2 .+ delta_k.^2))
    b = sum(dphi0_k ./ (sigma^2 .+ delta_k.^2))
    c = -0.5 * sum(dphi0_k.^2 ./ (sigma^2 .+ delta_k.^2))
    d = -0.5 * sum(log1p.((sigma ./ delta_k).^2))

    return a, b, c, d
end

function _evaluate_log_delta_phi_sigma_integrand(sigma::Float64, delta_phi::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; _internal_log_normalization_value::Union{Nothing, Float64} = nothing)
    # This function evaluates the integrand of the distribution p(δφ_n | data), at a given point σ, given the value of δφ_n (delta_phi)
    # _internal_log_normalization_value is a constant value that I will subtract from the log of the integrand, to regularize the pdf for numerical evaluation
    # In practice this amounts to dividing the integrand (and so the integral ~ pdf) by exp(_internal_log_normalization_value), which is a constant prefactor
    # This is perfectly allowed, since the posterior distribution is not normalized anyway

    if sigma < 0
        return 0.0
    end

    a, b, c, d = _evaluate_delta_phi_quantities_abcd(sigma, dphi0_k, delta_k)

    # println("a = ", a, ", b = ", b, ", c = ", c, ", d = ", d)

    log_integrand = -0.5 * (a * delta_phi^2 - 2. * b * delta_phi - (b*sigma)^2)/(1 + a * sigma^2) + c + d - 0.5 * log1p(a * sigma^2)

    if(_internal_log_normalization_value !== nothing)
        log_integrand -= _internal_log_normalization_value
    end

    # println("log_integrand (in sigma = ", sigma , ")= ", log_integrand)

    return log_integrand
end

function _evaluate_delta_phi_pdf(delta_phi::Float64, dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; _internal_log_normalization_value::Union{Nothing, Float64} = nothing)
    # This function estimates the integral over sigma, from 0 to +infinity, of the distribution p(δφ_n | data), given the value of δφ_n (delta_phi)
    # I evaluate the integral using the Gauss-Kronrod quadrature, variables the QuadGK package

    # Since this function is somewhat difficult to integrate, for improved accuracy I split the integral in three regions, using a guess for the spread of the distribution as a gauge of the order of magnitude over which the integrand varies
    guess_distributuion_spread = mean(delta_k)

    result1, err1 = quadgk(sigma -> exp(_evaluate_log_delta_phi_sigma_integrand(sigma, delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)), 0., guess_distributuion_spread, rtol = 1e-5, atol = 0.)
    result2, err2 = quadgk(sigma -> exp(_evaluate_log_delta_phi_sigma_integrand(sigma, delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)), guess_distributuion_spread, 5. *guess_distributuion_spread, rtol = 1e-5, atol = 0.)
    result3, err3 = quadgk(sigma -> exp(_evaluate_log_delta_phi_sigma_integrand(sigma, delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)), 5. *guess_distributuion_spread, + Inf, rtol = 1e-5, atol = 0.)

    result = result1 + result2 + result3
    err = sqrt(err1^2 + err2^2 + err3^2)

    # println("Integral result, for delta_phi = ", delta_phi , ": (", result1, " ± ", err1 , ") + (",  result2, " ± ", err2, ") + (", result3, " ± ", err3, ")")
    # println("Final integral result, for delta_phi = ", delta_phi , ": ", result, " ± ", err, ", relative error: ", err / abs(result))

    return result
end

# # Old function definition, with "agnostic" prior for delta_phi: not used since the convergence is really slow (could try to make the interval somewhat smaller, but probably will not be enough)
# function obtain_samples_delta_phi_pdf(dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; uniform_prior_range::NTuple{2, Union{Nothing, Float64}} = (nothing, nothing), n_samples::Int64 = 5000, burn_in::Int64 = 1000)
#     # This function samples from the distribution p(δφ_n | data) for a given set of parameters φ_k , Δ_k

#     # Evaluate a constant normalization factor, needed to regularize the pdf for numerical evaluation (I will subtract this -constant- value from the log of the integrand of pdf, so in practice I am dividing the pdf by exp(_internal_log_normalization_value), a constant prefactor)
#     # This is allowed, since the posterior distribution is not normalized anyway

#     # I estimate the order of magnitude of the spread of the integrand (2D function in sigma and delta_phi) by evaluating the mean of delta_k and abs(dphi0_k)
#     println("mean(delta_k) = ", mean(delta_k), ", mean(abs.(dphi0_k)) = ", mean(abs.(dphi0_k)))
#     # I create a 2x2 grid of 110 points, uniformly distributed in log space, for sigma and delta_phi
#     sigma_proposed_values = vcat(0., 10. .^ range(log10(1e-7*mean(delta_k)), log10(1e5*mean(delta_k)), 9))
#     delta_phi_proposed_values = vcat( (-1. * (10. .^ range(log10(1e-7*mean(abs.(dphi0_k))), log10(1e5*mean(abs.(dphi0_k))), 5))), 10. .^ range(log10(1e-7*mean(abs.(dphi0_k))), log10(1e5*mean(abs.(dphi0_k))), 5), 0.)

#     _internal_log_normalization_value = -Inf

#     # I evaluate the log of the integrand at the proposed points, to find the maximum of the log of the integrand
#     # which will be a constant that I will use to normalize the pdf, to improve the numerical stability of the integration
#     for sigma in sigma_proposed_values
#         for delta_phi in delta_phi_proposed_values
#             _internal_log_normalization_value_temp = _evaluate_log_delta_phi_sigma_integrand(sigma, delta_phi, dphi0_k, delta_k; _internal_log_normalization_value = nothing)
#             if _internal_log_normalization_value_temp > _internal_log_normalization_value
#                 _internal_log_normalization_value = _internal_log_normalization_value_temp
#                 # println("New internal log normalization value: ", _internal_log_normalization_value, " for sigma = ", sigma, ", delta_phi = ", delta_phi)
#             end    
#         end
#     end

#     println("Chosen internal log normalization value: ", _internal_log_normalization_value)


#     # _internal_log_normalization_value = _evaluate_log_delta_phi_sigma_integrand(0.0, initial_delta_phi_sampling_point, dphi0_k, delta_k; _internal_log_normalization_value = nothing)
#     # println("Internal log normalization value: ", _internal_log_normalization_value)
#     # _internal_log_normalization_value = _evaluate_log_delta_phi_sigma_integrand(mean(delta_k), initial_delta_phi_sampling_point, dphi0_k, delta_k; _internal_log_normalization_value = nothing)
#     # println("Internal log normalization value: ", _internal_log_normalization_value)
#     #     _internal_log_normalization_value = _evaluate_log_delta_phi_sigma_integrand(mean(delta_k), mean(dphi0_k), dphi0_k, delta_k; _internal_log_normalization_value = nothing)
#     # println("Internal log normalization value: ", _internal_log_normalization_value)
#     #_internal_log_normalization_value= 0.
    
#     # Now I will try to estimate a range for the uniform prior for delta_phi, based on the data
#     # In practice I will sample the function at arbitrary points, and chose the uniform range in the region where the returned integral is identically zero
#     # I assume no bimodality in the distribution or strange behavior...
#     delta_phi_proposed_values_prior = vcat(reverse(- 1. * (10. .^ range(log10(1e-5*mean(abs.(dphi0_k))), log10(1e10*mean(abs.(dphi0_k))), 15))), 10. .^ range(log10(1e-5*mean(abs.(dphi0_k))), log10(1e10*mean(abs.(dphi0_k))), 15))
#     probability_values = [_evaluate_delta_phi_pdf(x, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value) for x in delta_phi_proposed_values_prior]
#     # I will find the first and the last index where the probability is non-zero, and use that to define the range of the uniform prior
#     println(hcat(delta_phi_proposed_values_prior, probability_values))
#     first_nonzero_index = findfirst(probability_values .> 0.0)
#     last_nonzero_index = findlast(probability_values .> 0.0)
#     if first_nonzero_index === nothing || last_nonzero_index === nothing
#         error("No non-zero probability values found in the proposed range for delta_phi. Try to broaden the search range, or refine the grid.")
#     end
#     # I take just the maximum value of either, for symmetry, multiply by 10, for safety, and round them to the largest power of 10, for simplicity
#     uniform_prior_limit = 10. ^ (ceil(log10(10. * maximum(abs.((delta_phi_proposed_values_prior[first_nonzero_index], delta_phi_proposed_values_prior[last_nonzero_index]))))))

#     uniform_prior_range = (-uniform_prior_limit, uniform_prior_limit)
#     println("Uniform prior range for delta_phi (automatically chosen): ", uniform_prior_range)

#     @model function sample_model(dphi0_k, delta_k)
#         # Using a broad uniform prior for delta_phi, so mathematically I am still sampling from the posterior distribution p(δφ_n | data), if the prior is broad enough
#         delta_phi ~ Uniform(uniform_prior_range...)
        
#         # Add the posterior as a likelihood factor
#         pdf_value = _evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)
        
#         # Add log(pdf_value) to the total log probability
#         Turing.@addlogprob! log(pdf_value)
#     end

#     model = sample_model(dphi0_k, delta_k)

#     #Use the Metropolis-Hastings algorithm (or actually the no U-turns sampler, but I did not implement automatic differentiation here yet) to sample from the posterior distribution
#     chain = sample(model, MH(), n_samples + burn_in) #, init_params = (delta_phi = initial_delta_phi_sampling_point,))

#     #Remove the burn-in samples
#     samples = chain[:delta_phi][burn_in+1:end]

#     # Plot the chain and save the plot for debugging
#     plot(samples, title = "Samples from the posterior distribution p(δφ_n | data)", xlabel = "Sample index", ylabel = "δφ_n value")
#     savefig("samples_posterior_distribution.png")

#     return samples
# end



# New definition of the function, which uses information from the analytical form of the conditioned pdf (with sigma = 0)

"""
Function that evaluates the posterior distribution p(δφ_n | data) for the beyond GR deformation coefficient δφ_n, given the parameters φ_k , Δ_k, for k in (1,...,N) measurements of individual GW events.

Since this is an unnormalized 1D pdf (in delta_phi), I perform an MCMC to sample from the distribution. Then it will be possible to produce a violin plot with these samples, to plot the posterior distribution. 
Care should be taken of removing the burn-in samples and making sure that the chain has converged. 
I will produce trace plots and histograms of the samples, to check the convergence and the distribution of the samples.
"""
function obtain_samples_delta_phi_pdf(dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; n_samples::Int64 = 5000, burn_in::Int64 = 1000, debug_folder_name::Union{Nothing, String} = nothing, pnorder::Union{Nothing, String} = nothing)
    # This function samples from the distribution p(δφ_n | data) for a given set of parameters φ_k , Δ_k
    
    
    # Settings
    # Maximum number of MCMC runs to perform, to avoid infinite loops in case of problems with the prior range
    max_MCMC_runs = 4 
    # I will check that no points of the chain are near the boundaries of the prior range (as set by prior_fraction_to_be_empty_per_side, see also below). Otherwise, I issue a warning.
    prior_fraction_to_be_empty_per_side = 1 / 5.    


    # I estimate the order of magnitude of the maximum and the spread of the pdf (which I assume to be well behaved, and not multi-modal)
    # by using the analytical form of the pdf for the conditioned (sigma = 0) case (which is a normal distribution).
    # I will instaed estimate the range for the bulk of the integral over sigma by the mean of delta_k
    mu_conditioned_case, std_conditioned_case = getMuDistTIGER(dphi0_k, delta_k)
    println("mu_conditioned_case = ", mu_conditioned_case , ", std_conditioned_case = ", std_conditioned_case, " , mean(delta_k) = ", mean(delta_k), ", mean(abs.(dphi0_k)) = ", mean(abs.(dphi0_k)))

    # To increase the numerical accuracy and stability, I will also evaluate a constant normalization factor, needed to regularize the pdf for numerical evaluation 
    # (I will subtract this -constant- value from the log of the integrand of pdf, so in practice I am dividing the pdf by exp(_internal_log_normalization_value), a constant prefactor)
    # This is perfectly allowed, since the posterior distribution is not normalized anyway
    # Previously I evaluated this quantity over a grid of several points (in the 2D space in sigma and delta_phi), but now I will just evaluate it at a single point, informed by the conditioned distribution - hopefully this should suffice
    _internal_log_normalization_value = _evaluate_log_delta_phi_sigma_integrand(0., mu_conditioned_case, dphi0_k, delta_k; _internal_log_normalization_value = nothing)

    println("Internal log normalization value: ", _internal_log_normalization_value)


    # Since I cannot use a more efficient sampler (like Hamiltonian MonteCarlo or NUTS, since quadgk seems to break automatic differentiation), I will at least cache the last two computed value of _evaluate_delta_phi_pdf.
    # I do this because it seems that the MH() sampler in Turing does not cache the last computed values of the pdf, and so it recomputes them even when the last proposal was discarded, and so no move is made.
    # Since we are using the MH algorithm, I expect it will be enough to cache just the last two computed values of the pdf.
    local local_cache_evaluate_delta_phi_pdf = Tuple{Float64, Vector{Float64}, Vector{Float64}, Union{Nothing, Float64}, Float64}[]
    function cached_evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k; _internal_log_normalization_value=nothing)
        pdf_value = nothing
        if( length(local_cache_evaluate_delta_phi_pdf) > 1 )
            # I have cached enough values, I can start looking up the cached values
            # I assume to have to compare only the oldest value out of two. I also pop the first element, to keep the cache size small (two elements at most)
            values_two_evaluations_ago = popfirst!(local_cache_evaluate_delta_phi_pdf)
            if values_two_evaluations_ago[1:4] == (delta_phi, dphi0_k, delta_k, _internal_log_normalization_value)
                # The last cached value is the one I need
                pdf_value = values_two_evaluations_ago[5]
            else
                # I will check also the last element in the cache, to see if it matches
                if local_cache_evaluate_delta_phi_pdf[end][1:4] == (delta_phi, dphi0_k, delta_k, _internal_log_normalization_value)
                    # The last cached value is the one I need
                    pdf_value = local_cache_evaluate_delta_phi_pdf[end][5]
                else
                    # I have to compute the pdf value
                    pdf_value = _evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k; _internal_log_normalization_value=_internal_log_normalization_value)
                end
            end
        else
            pdf_value = _evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k; _internal_log_normalization_value=_internal_log_normalization_value)            
        end

        push!(local_cache_evaluate_delta_phi_pdf, (delta_phi, dphi0_k, delta_k, _internal_log_normalization_value, pdf_value))

        return pdf_value
    end

    # Actually, to further automatize the MCMC, I will performa while loop, where I widen the prior range, if some samples fell to close to the prior boundaries during the previous MCMC sampling 
    MCMC_run_index = 0
    samples = nothing
    uniform_prior_range = nothing

    while true 
        MCMC_run_index += 1
        if MCMC_run_index > max_MCMC_runs
            @warn "Performed already ($max_MCMC_runs) different runs of the MCMC, but some samples are still too close to the boundaries the prior range (within prior_fraction_to_be_empty_per_side = $(prior_fraction_to_be_empty_per_side) from the border)! This may indicate that the prior range is too small, or that the 'marginalized' posterior distribution is not well approximated by the conditioned pdf. You should increase uniform_prior_multiplier, to make sure that the prior range is large enough to encompass the bulk of the posterior distribution, and so not to bias the sampling process. Stopping the MCMC loop in any case."
    
        end
        println("MCMC run number: ", MCMC_run_index, " (out of a maximum of ", max_MCMC_runs, ")")

        # Now I will determine the range of the prior distribution, necessary for the MCMC sampling
        # If it is too wide, the converge will be too slow, so I try to pick a reasonable prior range, and if it is too small, it will be made larger in the next iteration
        if MCMC_run_index == 1
            # Now I will try to estimate a range for the uniform prior for delta_phi. 
            # Being the first run, I have no information on the 'marginalized' pdf itself, so I base my guess for the prior range on the previous displayed information, about the given quantities and the 'conditioned' pdf
            # I assume that the conditioned and this "marginalized" pdf are quite similar
            # Therefore for the prior I choose a uniform interval, centered around mu_conditioned_case, and with width 2 * (uniform_prior_multiplier * std_conditioned_case)
            #uniform_prior_multiplier = 20.
            uniform_prior_multiplier = 30.
            # I will also multiply the prior range by the factor 1. / (1. - 2. * prior_fraction_to_be_empty_per_side ), so that I take into account the 'unusable' prior range
            uniform_prior_multiplier = uniform_prior_multiplier * (1. / (1. - 2. * prior_fraction_to_be_empty_per_side ))

            uniform_prior_range = (mu_conditioned_case - (uniform_prior_multiplier * std_conditioned_case), mu_conditioned_case + (uniform_prior_multiplier * std_conditioned_case))
            # I expect the spread of this "marginalized" pdf to be larger than the conditioned one, but by considering a factor of 20 [30] (uniform_prior_multiplier) I should be safe 
            # (In fact I will still be able to correctly sample points within about the 5 sigma interval, even if the spread of this "marginalized" pdf is 4 times larger than the conditioned one)
            # For safety this value should be made larger, but this decreases the convergence speed of the MCMC sampling
            # Also the mean of the distribution may be slightly different
            # I could also perform some checks right now, for example estimating the maximum of the pdf in the prior range, and then checking that at the boundaries of the prior range the probability density is much lower than the maximum
            # However I will directly sample via MCMC, and instead look at the resulting chain
            # => Therefore, a posteriori, as a rule of thumb, one should always chech that recovered samples are always well within the bounds of the prior! I will do so below, producing plots and performing an automatic check.

        elseif MCMC_run_index > 1
            # If this is not the first run, it means some samples were too close to the boundaries of the prior range (as set by prior_fraction_to_be_empty_per_side, see also above)
            # In this case then, I will widen the prior range, and repeat the MCMC sampling
            
            # I will find the maximum value of the samples, and also the mean and the sigma
            max_sample_value = maximum(samples)
            min_sample_value = minimum(samples)
            mean_sample_value = mean(samples)
            std_sample_value = std(samples)
            
            # I now set the prior range to be 5 * std_sample_value * 2^(MCMC_run_index), centered around mean_sample_value
            uniform_prior_multiplier = 5. * 2. ^ (MCMC_run_index) * (1. / (1. - 2. * prior_fraction_to_be_empty_per_side ))
            # but if min_sample_value or max_sample_value are already too close to the boundaries, I will continue to the next cycle, to further widen the prior range
            uniform_prior_range = (mean_sample_value - (uniform_prior_multiplier * std_sample_value), mean_sample_value + (uniform_prior_multiplier * std_sample_value))

            if (min_sample_value < (mean_sample_value - (uniform_prior_multiplier * std_sample_value * (1. - 2. * prior_fraction_to_be_empty_per_side )))) || (max_sample_value > (mean_sample_value + (uniform_prior_multiplier * std_sample_value * (1. - 2. * prior_fraction_to_be_empty_per_side ))))
                println("Samples from the previous MCMC run would alreby be too close to the boundaries of the new prior range (uniform_prior_range = $uniform_prior_range): skipping this MCMC run altogether, and further widening the prior range.")
                continue # Go to the next iteration of the while loop, to widen the prior range
            end
        else
            @error "MCMC run index is not valid: ($MCMC_run_index). Aborting."
            return nothing
        end
        
        println("Uniform prior range for delta_phi (automatically chosen): ", uniform_prior_range)
        



        # Perform the MCMC sampling from the posterior distribution p(δφ_n | data)

        @model function sample_model(dphi0_k, delta_k)
            # Using a broad uniform prior for delta_phi, so mathematically I am still sampling from the posterior distribution p(δφ_n | data), if the prior is broad enough
            delta_phi ~ Uniform(uniform_prior_range...)
            
            # Add the posterior as a likelihood factor
            # pdf_value = _evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)
            # I use the cached version of the pdf evaluation, to speed up the MCMC sampling
            pdf_value = cached_evaluate_delta_phi_pdf(delta_phi, dphi0_k, delta_k, _internal_log_normalization_value = _internal_log_normalization_value)
            
            # Add log(pdf_value) to the total log probability
            Turing.@addlogprob! log(pdf_value)
        end

        model = sample_model(dphi0_k, delta_k)

        #Use the Metropolis-Hastings algorithm (I would like to use the no U-turns sampler, but I did not implement automatic differentiation here yet) to sample from the posterior distribution
        chain = sample(model, MH(), n_samples + burn_in) #, init_params = (delta_phi = initial_delta_phi_sampling_point,))

        # For debugging purposes, I will plot the trace plot of the chain, and histogram the points
        if debug_folder_name !== nothing && pnorder !== nothing
            # Thread locking to avoid multiple threads trying to plot at the same time (plot backend is not thread-safe)
            lock(plot_lock) do
                # Create the folder if it does not exist
                mkpath(debug_folder_name)
                # For debugging purposes, I will plot the chain, with a vertical line where the burn in ends
                # Also, in the same plot, in the upper part of the figure, the samples from the posterior distribution p(δφ_n | data), with the uniform distribution range shown, and in the lower part the chain
                # By visual inspection, the chain should be stable after the burn-in period, and the samples should be well distributed within the prior usable_prior_range
                # Create the trace plot (chain)
                trace_plot = plot(chain[:delta_phi], title = "Trace plot for the posterior distribution p(δφ_n | data)", xlabel = "Sample index", ylabel = "δφ_n value", label = "Chain samples")
                vline!(trace_plot, [burn_in], label = "Burn-in end", color = :red, linestyle = :dash, lw = 3)

                # Create the histogram
                hist_xmin = uniform_prior_range[1] - 0.03 * (uniform_prior_range[2] - uniform_prior_range[1])
                hist_xmax = uniform_prior_range[2] + 0.03 * (uniform_prior_range[2] - uniform_prior_range[1])
                
                hist_plot = histogram(chain[:delta_phi][burn_in+1:end], title = "Posterior distribution p(δφ_n | data)", xlabel = "δφ_n value", ylabel = "Prob. density", label = "Posterior samples histogram (burn in removed)", legend = :topright, normalize = true, xlims = (hist_xmin, hist_xmax))
                # Shade outside of prior range
                vspan!(hist_plot, [hist_xmin, uniform_prior_range[1]], color = RGBA(0., 0., 0., 0.8), label = "Outside prior range")
                vspan!(hist_plot, [uniform_prior_range[2], hist_xmax], color = RGBA(0., 0., 0., 0.8), label = "")
                vspan!(hist_plot, [uniform_prior_range[1], uniform_prior_range[1] + prior_fraction_to_be_empty_per_side * (uniform_prior_range[2] - uniform_prior_range[1])], color = RGBA(0.7, 0.8, 1.0, 0.4), label = "Inside prior range, but wanted to be empty")
                vspan!(hist_plot, [uniform_prior_range[2] - prior_fraction_to_be_empty_per_side * (uniform_prior_range[2] - uniform_prior_range[1]), uniform_prior_range[2]], color = RGBA(0.7, 0.8, 1.0, 0.4), label = "")
                
                # Overlay the distribution of the conditioned case
                plot!(hist_plot, Normal(mu_conditioned_case, std_conditioned_case), color=:yellow, lw=4, label="Conditioned (sigma = 0) pdf for reference")

                # Draw vertical lines for the prior range
                vline!(hist_plot, [uniform_prior_range[1], uniform_prior_range[2]], label = "Uniform prior range", color = :orange, linestyle = :dash, lw = 3)
                
                # Combine both plots in a vertical layout (2 rows, 1 column)
                plt = plot(hist_plot, trace_plot, layout = @layout([a; b]), size = (2000, 2000))

                # Save the combined figure
                savefig(plt, joinpath(debug_folder_name, "MCMC_posterior_distribution_pn_$(pnorder)_run$(MCMC_run_index).pdf"))
                println("Saved the trace plot of the chain and histogram of the posterior samples to: ", joinpath(debug_folder_name, "MCMC_posterior_distribution_pn_$(pnorder)_run$(MCMC_run_index).png"))      
            end
        else
            println("No debug folder and pnorder provided, skipping the debug plot of the chain.")
        end


        #Remove the burn-in samples
        samples = chain[:delta_phi][burn_in+1:end]

        # Check that no points of the chain are near the boundaries of the prior range (I would like no points at all to fall within the first 1/5 or within last 1/5 of the prior range, set by prior_fraction_to_be_empty_per_side before)
        # If this is close enough to a gaussian distribution (assuming the mean to be centered in the interval...), and you are sampling for example N (e.g. 10000) points, 
        # then the probability of having by chance at least one point within the first 1/5 or last 1/5 of the prior range is given by p_outside_total, which can becomputed here below:
        #
        # using Distributions
        # N_samples = 10000
        # sigma_value = 1.0         # sigma of the distribution (arbitrary units)
        # width_prior = 2. * 20.    # 2 * uniform_prior_multiplier (in the same arbitrary units of sigma_value)
        # prior_fraction_to_be_empty_per_side = 1 / 5.
        # usable_prior_range = (width_prior * (1. - 2. * prior_fraction_to_be_empty_per_side) ) / 2.
        # # Probability that a single point lies within delta_phi lies within ± usable_prior_range
        # p_within = cdf(Normal(0, sigma_value), usable_prior_range) - cdf(Normal(0, sigma_value), - usable_prior_range)
        # # Probability that at least one point  delta_phi lies outside ± usable_prior_range
        # p_outside_total = 1. - (p_within)^N_samples
        #
        # As an example, for N_samples = 10000, with uniform_prior_multiplier = 20, the spread of the conditioned distribution could be 2 times larger than the one of the conditioned distribution, 
        # and still you would have about 0.001% chance of having at least one point within the first 1/5 or last 1/5 of the prior range. If instead the mean is translated, or the sigma is larger than about twice,
        # it becomes highly likely that you obtain at least one point within the first 1/5 or last 1/5 of the prior range, and so a warning here below.
        # In that case nonetheless it would be wise to increase the uniform_prior_multiplier, to make sure that the prior range is large enough to encompass the bulk of the posterior distribution, and so not to bias the sampling process.

        if any(samples .< uniform_prior_range[1] + prior_fraction_to_be_empty_per_side * (uniform_prior_range[2] - uniform_prior_range[1])) || any(samples .> uniform_prior_range[2] - prior_fraction_to_be_empty_per_side * (uniform_prior_range[2] - uniform_prior_range[1]))
            println("Some MCMC samples are close to the boundaries the prior range (within prior_fraction_to_be_empty_per_side = $(prior_fraction_to_be_empty_per_side) from the border)! This may indicate that the prior range is too small, or that the 'marginalized' posterior distribution is not well approximated by the conditioned pdf. I will repeat the MCMC sampling, widening the prior range.")
        else
            println("All MCMC samples obtained are further than required from the boundaries of the prior range (prior_fraction_to_be_empty_per_side = $(prior_fraction_to_be_empty_per_side)). Returning the samples.")
            break
        end
        
    end

    return collect(samples)
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

function obtain_samples_delta_phi_pdf_conditioned(dphi0_k::Vector{Float64}, delta_k::Vector{Float64}; n_samples::Int64 = 5000)
    # This function samples from the distribution p(δφ_n | data) for a given set of parameters φ_k , Δ_k, with the posterior distribution conditioned on σ=0
    # I know the analytical form of the distribution, which is a Normal distribution: then I will just sample from it
    
    meanMuDist, stdMuDist = getMuDistTIGER(dphi0_k, delta_k)
    dist = Normal(meanMuDist, stdMuDist)      # Create a Normal distribution
    samples = rand(dist, n_samples)

    return samples
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

function obtain_conditioned_upper_bounds(n_events::Integer, global_index_network::Vector{Bool}, dphi0_k_full::Vector{Float64}, delta_k_full::Vector{Float64}; averageOverSeveralRealizations::Bool = true, numberOfEventsSingleRealization::Union{Nothing, Int64} = nothing, n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events::Bool = false, printEventsAsHorizontalLinesOrDensityPlot::Bool = true, use_all_n_events_for_single_event_sample_distribution::Bool = true)
    if n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events
        println("Considering directly observed events (i.e. not performing a random draw (which would introduce 'poissonian' noise - actually binomially distributed) from the catalog)! Therefore, selecting directly n_events = $(n_events) from the observed ones! This may not be what you want!")
        # If the number of events refers to the observed events, then we select the events to be used (in practice dphi0_k and delta_k) directly only from observed events (i.e. with global_index_network == 1)
        
        # extract only valid data-points via global index
        if n_events > sum(global_index_network)
            println("Warning: There are at most $(sum(global_index_network))-observed events available")
            n_events_used = sum(global_index_network)
        else
            n_events_used = n_events
            println("Observed vents used: $(n_events_used)")
        end

        # Not meaningful here, since it will be just an array of ones, but I keep it for consistency because in the other if branch it contains non trivial information
        global_index_network_effective = global_index_network[global_index_network]

        dphi0_k = dphi0_k_full[global_index_network]
        delta_k = delta_k_full[global_index_network]
    else
        println("Drawing events from the catalog (so with associated poissonian noise  - actually binomially distributed - on the counting of observed events)! So drawing a number of observed events starting from n_events = $(n_events) obtained from the catalog!")
        if n_events > length(global_index_network)
            println("Warning: There are at most $(length(global_index_network))-events in the catalog available")
            n_events_used = length(global_index_network)
        else
            n_events_used = n_events
            println("Events used (from the catalog): $(n_events_used)")
        end
        
        # In this case it is useful, since it contains the information about which events will be considered as 'observed', after the draw from the catalog, and the several realizations
        global_index_network_effective = global_index_network

        dphi0_k = dphi0_k_full
        delta_k = delta_k_full
    end

    global_index_network_effective = global_index_network_effective[1:n_events_used]
    delta_k = delta_k[1:n_events_used]
    dphi0_k = dphi0_k[1:n_events_used] 

    if averageOverSeveralRealizations
        if numberOfEventsSingleRealization === nothing
           throw(ArgumentError("numberOfEventsSingleRealization must be specified when averageOverSeveralRealizations is true."))
        end
        numberOfRealization = Int(floor(n_events_used / numberOfEventsSingleRealization))
        if numberOfRealization == 0
            @error "Not enough events to perform the average over several realizations. Please decrease the numberOfEventsSingleRealization in the config file, or increase the number of total events."
        end
        n_events_used = numberOfEventsSingleRealization * numberOfRealization

        println("Splitting the dataset into $(numberOfRealization) realizations of $(numberOfEventsSingleRealization) events")
        println("Number of events used effectively: $(n_events_used)")
        
    else
        numberOfRealization = 1
        numberOfEventsSingleRealization = n_events_used
    end

    global_index_network_effective = reshape(global_index_network_effective[1:n_events_used], numberOfRealization,  numberOfEventsSingleRealization)
    dphi0_k = reshape(dphi0_k[1:n_events_used], numberOfRealization,  numberOfEventsSingleRealization)
    delta_k = reshape(delta_k[1:n_events_used], numberOfRealization, numberOfEventsSingleRealization)

    vectorUpperLimits = zeros(numberOfRealization)
    dphi0_k_realization = Vector{Vector{Float64}}(undef, numberOfRealization)
    delta_k_realization = Vector{Vector{Float64}}(undef, numberOfRealization)
    number_events_single_realization = zeros(Int64, numberOfRealization)

    for realization_index in 1:numberOfRealization
        # I will select only the events which pass the selection cuts (SNR, invertible Fisher, eventually inspiral SNR thresholds), as indicated by global_index_network_effective
        dphi0_k_realization[realization_index] = dphi0_k[realization_index, global_index_network_effective[realization_index, :]]
        delta_k_realization[realization_index] = delta_k[realization_index, global_index_network_effective[realization_index, :]]
        number_events_single_realization[realization_index] = sum(global_index_network_effective[realization_index, :])
    end

    for realization_index in 1:numberOfRealization

        # Now I iterate over each realization
        if number_events_single_realization[realization_index] == 0
            # No events in realization, so I will skip this realization and set the upper limit to NaN
            vectorUpperLimits[realization_index] = NaN
            continue
        end

        #Evaluate the 90% upper limit given this single realization
        vectorUpperLimits[realization_index] = get90PctUpperLimitMuDistTIGER(dphi0_k_realization[realization_index], delta_k_realization[realization_index], 0.9)

    end
    

    upperLimitSingleEventsTemp = nothing # This will be used to store the upper limits for the single events, if requested

    #Save the single events 90% upper limits (eventually for a given realization), if requested
    if printEventsAsHorizontalLinesOrDensityPlot
        # The hierarchical analysis, conditioned on sigma = 0, implemented in getMuDistTIGER, should work just fine even if working with a single event

        if use_all_n_events_for_single_event_sample_distribution
            # I save all the single events upper limits for the first realization, so that if plotted as a density plot the fluctuactions are smaller
            #I flatten dphi0_k and delta_k over the realization_index with vec, so that I can use all the event
            upperLimitSingleEventsTemp = map((x, y) -> get90PctUpperLimitMuDistTIGER([x], [y], 0.9), vcat(dphi0_k_realization...), vcat(delta_k_realization...))
        else
            # I use only a single realization to plot the single events upper limits
            # I look for the first realization_index which has at least one event
            first_realization_with_events = findfirst(x -> x > 0, number_events_single_realization)
            upperLimitSingleEventsTemp = map((x, y) -> get90PctUpperLimitMuDistTIGER([x], [y], 0.9), dphi0_k_realization[first_realization_with_events, :], delta_k_realization[first_realization_with_events, :])
        end
        
    end

    return vectorUpperLimits, upperLimitSingleEventsTemp, number_events_single_realization, numberOfEventsSingleRealization
end
