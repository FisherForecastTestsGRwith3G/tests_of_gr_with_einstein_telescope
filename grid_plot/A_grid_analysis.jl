using HDF5
using TOML
#using Plots
using Trapz
#using LaTeXStrings
using KernelDensity, Statistics
using Serialization
#using Turing
#using JSON
# include required scripts 
include("../hierachical_combination/hierarchical_distribution.jl")
include("_config_parser_grid.jl")
# include("_conditioned_plot_utils.jl")
include("_grid_utils.jl")

# obtain config file from input
# user_configs = getUserConfigs()
# config_file_name = "config_files/"*user_configs["default_config"]

# check if we should run the simulation
needToRun = true
# needToRun = false
# if length(ARGS) > 0
#     if ARGS[1] == "1"
#         needToRun = true
#         println("The simulation will be run!")
#     elseif ARGS[1] == "0"
#         needToRun = false
#         println("The simulation will not be run!")
#     else
#         @warn "There was a first input argument handed. But it was $(ARGS[1]) and not \"1\" or \"0\" and thus ignored"
#     end
# end

# if length(ARGS) > 1
#     config_file_name = ARGS[2]
# else
#     println("No config file handed, using default!")
# end
config_file_name = ARGS[1]
println("Using config file: $(config_file_name)\n")

configs = read_config_grid(config_file_name)


mu_vec = configs["mu_vec"]
sigma_vec = configs["sigma_vec"]
PN_string = configs["pn_orders"][1]
network = configs["network"]
n_events = configs["n_events"]
gridSize = length(mu_vec) * length(sigma_vec) 

if needToRun
    # run the simulation
    println("Running the simulation...")
    local idx = 0
    # local skips = configs["skips"]

    for i in 1:length(mu_vec)
        for j in 1:length(sigma_vec)
            idx += 1
            println("Running simulation for index $(idx) out of $(gridSize)")
            # skip the (0.,0.) case
            if mu_vec[i] == 0. && sigma_vec[j] == 0.
                continue
            end

            # if skips > 0
            #     skips -=1
            #     continue
            # end
            # modify the config file
            header = "grid_PN_$(PN_string)_n_$(idx)_mu_" * string(mu_vec[i]) * "_sigma_" * string(sigma_vec[j])

            mu = mu_vec[i]
            sigma = sigma_vec[j]
            mkpath("config_files/grid/")
            config_file_name_out = "config_files/grid/config_A_"*header*".toml"

            modify_configs(config_file_name, config_file_name_out, header, mu, sigma, PN_string, network, n_events)
            # run the simulation
            
            run(`julia --project=. scripts_to_run/A_fisher_analysis.jl $(config_file_name_out)`)
            println("Simulation for index $(idx) out of $(gridSize) done!")

        end
    end
else
    println("The simulation will not be run!")
end
