using Random
using Base.Threads
using Statistics
using TOML
using Interpolations

function modify_configs(config_file_name, config_file_name_out, header, mu, sigma, PN, network, n_events)

    # Read the TOML file into a Julia dictionary
    config_dic = TOML.parsefile(config_file_name)

    delete!(config_dic["deviations"], "mu_vec")
    delete!(config_dic["deviations"], "sigma_vec")
    config_dev =config_dic["deviations"]
    config_cat = config_dic["catalog"]

    # Modify the values
    config_dev["mu"] = mu
    config_dev["sigma"] = sigma
    config_cat["n_events"] = n_events

    #config_dic["header"] = header
    #config_glob["pn_waveforms"] = [PN]
    #config_glob["network_list"] = network


    # Write the modified dictionary back to the file
    open(config_file_name_out, "w") do io
        TOML.print(io, config_dic)
    end
end




function bisection_method(f, a, b; tol=1e-6, max_iter=100, debug=true, args_f=())
    """
    Bisection method to find the root of a function.

    Parameters:
    - f: Function for which the root is to be found.
    - a: Lower bound of the interval.
    - b: Upper bound of the interval.
    - tol: Tolerance for convergence (default: 1e-6).
    - max_iter: Maximum number of iterations (default: 100).

    Returns:
    - root: Approximation of the root.
    - converged: Boolean indicating if the method converged.
    """       
    try 

        fa = f(a, args_f...)
        fb = f(b, args_f...)
        # if res > 0 it means that (0., 0.) is outside the 3 sigma level

        if  fa * fb > 0
            if debug
                println("The function must have opposite signs at the endpoints a and b.")
            end
            if fa > 0.
            
                return a, false
            else
                return b, false
            end
        end

        for i in 1:max_iter
            c = (a + b) / 2  # Midpoint
            fc = f(c, args_f...)

            if abs(fc) < tol || abs(b - a) < tol
                return c, true, i
            end

            if f(a, args_f...) * fc < 0
                b = c  # Root is in [a, c]
            else
                a = c  # Root is in [c, b]
            end
        end

        if debug
            println("Maximum iterations reached without convergence.")
        end
        return (a + b) / 2, false
    catch e
        # if debug
        #     println("Error in bisection method: ", e)
        # end
        #@error "ERROR: " exception=(err, catch_backtrace())
        rethrow(e)
        return NaN, false
    end
end


function reshuffling_bisection(n_reshuffling, ff, dphi0_k, delta_k, mu, sigma)
    res = zeros(n_reshuffling)

    n_events = length(delta_k)
    @threads for i in 1:n_reshuffling

        p = randperm(n_events)

        dphi0_k_shuf = dphi0_k[p]
        delta_k_shuf = delta_k[p]

        res[i] = bisection_method(ff, 3, n_events, tol=1,  args_f=(dphi0_k_shuf, delta_k_shuf, sigma, mu))[1]
    end

    return median(res), res
end

    


function sigmas_from_center(pdf::Matrix{Float64}, x::Float64, y::Float64)
    # Find the center of the distribution
    nx, ny = size(pdf)
    cx, cy = nx ÷ 2, ny ÷ 2  # Assuming the center is at the middle of the grid

    # Calculate the standard deviation of the distribution
    mu_x = sum(pdf[i, j] * i for i in 1:nx, j in 1:ny) / sum(pdf)
    mu_y = sum(pdf[i, j] * j for i in 1:nx, j in 1:ny) / sum(pdf)
    std_x = sqrt(sum(pdf[i, j] * (i - mu_x)^2 for i in 1:nx, j in 1:ny) / sum(pdf))
    std_y = sqrt(sum(pdf[i, j] * (j - mu_y)^2 for i in 1:nx, j in 1:ny) / sum(pdf))

    # Calculate sigmas from the center
    sigmas_x = (x - mu_x) / std_x
    sigmas_y = (y - mu_y) / std_y

    return sigmas_x, sigmas_y
end


function wrapper_3sigma(n_events_used, dphi0_k, delta_k, center_mu, center_sig)

    n_events_used = Int(round(n_events_used))

    dphi0_k = dphi0_k[1:n_events_used] 
    delta_k = delta_k[1:n_events_used]

    ### calculate the hyper-parameter distribution
    # first estimate where to place it
    #center_mu = sum(dphi0_k) / n_events_used
    # center_sig = max(0, sqrt.(sum((center_mu .- dphi0_k).^2) ./  n_events_used))
    spread = sqrt.(1.0 ./ sum(1 ./ delta_k.^2) )

    k_spread = configs_B["k_spread"]
    mu_limit = (center_mu - k_spread*spread, center_mu + k_spread*spread)
    if !isnothing(configs_B["mu_lims"][pno]) # overwrite mu_lims
        mu_limit = (configs_B["mu_lims"][pno][1], configs_B["mu_lims"][pno][2])
    end 
    
    sig_limit = (max(0, center_sig - k_spread*spread), center_sig + k_spread*spread)
    if !isnothing(configs_B["sigma_lims"][pno]) # overwrite sig_lims
        sig_limit = (configs_B["sigma_lims"][pno][1], configs_B["sigma_lims"][pno][2])
    end 
    #sig_limit = (0., configs_B["sigma_lims"][pno][2])

    # println("mu_limit: ", mu_limit)
    # println("sig_limit: ", sig_limit)

    # calculate the ranges 
    mu_values = collect(LinRange(mu_limit[1], mu_limit[2], configs_B["n_points"]))
    sig_values = collect(LinRange(sig_limit[1], sig_limit[2], configs_B["n_points"]))
    p_mu_sig, p_sig, p_mu, n_tot, network_marg = hyperparamDistTIGER(mu_values, sig_values, dphi0_k, delta_k)

    # 2d interpolation
    itp = interpolate((mu_values, sig_values), p_mu_sig, Gridded(Linear()))
    p_mu_sig_interp = extrapolate(itp, 0.0)

    local level, p_GR
    try
        level = calcPercentileLvl( p_mu_sig, [0.9889])
        p_GR = p_mu_sig_interp(0., 0.)
    catch e
        println("Error in calculating percentiles: ", e)
        # If the percentile calculation fails, return 1, # indicating that the point (0., 0.) is outside the 3 sigma level.
        return 1.
    end

    # out of the 3 sigma interval if level > p_GR
    res = level[1] - p_GR
    # if res > 0 it means that (0., 0.) is outside the 3 sigma level
    return res

end

