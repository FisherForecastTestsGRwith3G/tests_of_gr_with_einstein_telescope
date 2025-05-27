using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# include required scripts 
include("_hierarchical_dist.jl")
include("_hierarchical_plotting_utils.jl")
include("_parse_config.jl")
include("_grid_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForZ(config_file_name)
configs_B, simulation_tag_B = readConfigForB(config_file_name)


# set up folder names
data_folder_name = user_configs["path_output"]*simulation_tag * "/"
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"

mu_vec = configs["mu"]
sigma_vec = configs["sigma"]
pno = configs["pn_waveforms"][1]
PN_name = pn_order_dic[pno][2]
network = configs["network_list"][1]
n_events = configs["n_events"]
gridSize = length(mu_vec) * length(sigma_vec) 

# get global index
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()


# process all the plots

idx = 0

for i in 1:length(mu_vec)
    for j in 1:length(sigma_vec)
        idx += 1
        println("Analyzing data for index $(idx) out of $(gridSize)")
        # skip the (0.,0.) case
        if mu_vec[i] == 0. && sigma_vec[j] == 0.
            continue
        end

        mu = mu_vec[i]
        sigma = sigma_vec[j]
        folder = data_folder_name * nn * "/pn_" * pn_order_dic[pno][2]
        filename = folder * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end
        file_name = data_grid_folder_name * "/global_indices.h5"
        h5open(file_name, "r") do gi_file
            global_index_fisher[network] = read(gi_file, "fisher")  
            global_index_snr[network] = read(gi_file, "snr")
            global_index_total[network] = read(gi_file, "total")
        end

        
        filename = data_grid_folder_name * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # extract only valid data-points via global index
        dphi0_k = data["dphi0_k"][global_index_total[network]]
        delta_k = data["delta_k"][global_index_total[network]]

        function wrapper_3sigma(n_events_used, dphi0_k, delta_k)

            dphi0_k = dphi0_k[1:n_events_used] 
            delta_k = delta_k[1:n_events_used]

            ### calculate the hyper-parameter distribution
            # first estimate where to place it
            center_mu = sum(dphi0_k) / n_events_used
            center_sig = max(0, sqrt.(sum(center_mu .- dphi0_k).^2 ./  n_events_used))
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

            # calculate the ranges 
            mu_values = collect(LinRange(mu_limit[1], mu_limit[2], configs["n_points"]))
            sig_values = collect(LinRange(sig_limit[1], sig_limit[2], configs["n_points"]))
            p_mu_sig, p_sig, p_mu, n_tot, network_marg = hyperparamDistTIGER(mu_values, sig_values, dphi0_k, delta_k)

            return p_mu_sig
    
            # calculate the 3 sigma upper limit
            # first estimate where to place it  

end
end

