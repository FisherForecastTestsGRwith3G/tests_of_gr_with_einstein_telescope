using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing
using JSON

# include required scripts 
include("_hierarchical_dist.jl")
include("_parse_config.jl")
include("_conditioned_plot_utils.jl")
include("_grid_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]

# check if we should run the simulation
needToRun = false
if length(ARGS) > 0
    if ARGS[1] == "1"
        needToRun = true
        println("The simulation will be run!")
    elseif ARGS[1] == "0"
        needToRun = false
        println("The simulation will not be run!")
    else
        @warn "There was a first input argument handed. But it was $(ARGS[1]) and not \"1\" or \"0\" and thus ignored"
    end
end

if length(ARGS) > 1
    config_file_name = ARGS[2]
else
    println("No config file handed, using default!")
end

println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForZ(config_file_name)

# set up folder names
data_folder_name = user_configs["path_output"]*simulation_tag*"/data/"
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"

mu_vec = configs["mu"]
sigma_vec = configs["sigma"]
PN = configs["pn_waveforms"][1]
PN_name = pn_order_dic[PN][2]
network_list = configs["network_list"]
n_events = configs["n_events"]
gridSize = length(mu_vec) * length(sigma_vec) 

if needToRun
    # run the simulation
    println("Running the simulation...")
    local idx = 0
    for i in 1:length(mu_vec)
        for j in 1:length(sigma_vec)
            idx += 1
            println("Running simulation for index $(idx) out of $(gridSize)")
            # skip the (0.,0.) case
            if mu_vec[i] == 0. && sigma_vec[j] == 0.
                continue
            end
            # modify the config file
            #header = "grid_PN_$(PN_name)_n_$(idx)_mu_" * string(mu_vec[i]) * "_sigma_" * string(sigma_vec[j])
            header = "grid_PN_$(PN_name)_n_$(idx)_mu_" * string(mu_vec[i]) * "_sigma_" * string(sigma_vec[j])

            mu = mu_vec[i]
            sigma = sigma_vec[j]
            config_file_name_out = "config_files/grid/config_A_"*header*".json"

            modify_configs(config_file_name, config_file_name_out, simulation_tag*"/"*header, mu, sigma, PN, network_list, n_events)
            # run the simulation
            run(`julia A_run_catalog_script.jl 1 $(config_file_name_out)`)
            println("Simulation for index $(idx) out of $(gridSize) done!")

        end
    end
else
    println("The simulation will not be run!")
end

# # get global index
# global_index_fisher = Dict()
# global_index_snr = Dict()
# global_index_total = Dict()
# for nn in configs["network_list"]
#     file_name = data_folder_name * nn *"/global_indices.h5"
#     h5open(file_name, "r") do gi_file
#         global_index_fisher[nn] = read(gi_file, "fisher")  
#         global_index_snr[nn] = read(gi_file, "snr")
#         global_index_total[nn] = read(gi_file, "total")
#     end
# end