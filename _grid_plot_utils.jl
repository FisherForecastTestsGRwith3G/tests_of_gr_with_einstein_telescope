


# function modify_configs(config_name, header, mu, sigma)
    
#     open(config_name, "w") do

#     # Modify the configs dictionary to set the values of mu and sigma
#     configs["mu"] = mu
#     configs["sigma"] = sigma
#     configs["header"] = header
#     # configs["simulation_tag"] = header * "_mu_" * string(mu) * "_sigma_" * string(sigma)

#     return configs
# end


# using JSON
function modify_configs(config_file_name, config_file_name_out, header, mu, sigma, PN, network_list, n_events)

    # Read the JSON file into a Julia dictionary
    config_dic = open(config_file_name,"r") do f
        config_dic = JSON.parse(f)
    end

    config_glob =config_dic["global"]
    config_a = config_dic["a_specific"]

    # Modify the values
    config_a["mu"] = mu
    config_a["sigma"] = sigma
    config_a["n_events"] = n_events

    config_dic["header"] = header
    config_glob["pn_waveforms"] = [PN]
    config_glob["network_list"] = network_list


    # Write the modified dictionary back to the file
    open(config_file_name_out, "w") do io
        JSON.print(io, config_dic, 4)
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
        if f(a, args_f...) * f(b, args_f...) > 0
            if debug
                println("The function must have opposite signs at the endpoints a and b.")
            end
            return a, false
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
        if debug
            println("Error in bisection method: ", e)
        end
        return NaN, false
    end
end

using Statistics

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

# # Example usage:
# # Assuming you have a 2D pdf matrix `pdf` and coordinates `x_point`, `y_point`
# # pdf = your_2d_pdf_matrix
# # x_point = x-coordinate of the point of interest
# # y_point = y-coordinate of the point of interest

# # Calculate sigmas from center for a specific point
# sigmas_x, sigmas_y = sigmas_from_center(pdf, x_point, y_point)
# println("Sigmas from center (x-direction): ", sigmas_x)
# println("Sigmas from center (y-direction): ", sigmas_y)
