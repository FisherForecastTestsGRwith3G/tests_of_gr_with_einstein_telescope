using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing
using Base.Threads


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
n_median = configs["n_median"]
gridSize = length(mu_vec) * length(sigma_vec) 

println("\n Anlayzing grid composed of $gridSize grid points\n")

# get global index
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()


# process all the plots


idx=0
res = zeros((length(mu_vec), length(sigma_vec)))
for i in 1:length(mu_vec)
    for j in 1:length(sigma_vec)
        global idx += 1
        println("\n Analyzing data for index $(idx) out of $(gridSize)")
        # skip the (0.,0.) case
        if mu_vec[i] == 0. && sigma_vec[j] == 0.
            continue
        end

        mu = mu_vec[i]
        sigma = sigma_vec[j]
        header = "grid_PN_$(PN_name)_n_$(idx)_mu_" * string(mu) * "_sigma_" * string(sigma)
        data_grid_folder_name = data_folder_name * header * "/data/" * network# * "/pn_" * PN_name 

        println(" Mu: ", mu, " Sigma: ", sigma)
        file_name = data_grid_folder_name * "/global_indices.h5"
        h5open(file_name, "r") do gi_file
            global_index_fisher[network] = read(gi_file, "fisher")  
            global_index_snr[network] = read(gi_file, "snr")
            global_index_total[network] = read(gi_file, "total")
        end

        data_grid_folder_name_pno = data_folder_name * header * "/data/" * network * "/pn_" * PN_name
        filename = data_grid_folder_name_pno * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # extract only valid data-points via global index
        dphi0_k = data["dphi0_k"][global_index_total[network]]
        delta_k = data["delta_k"][global_index_total[network]]

        res[i,j] = reshuffling_bisection(n_median, wrapper_3sigma, dphi0_k, delta_k, sigma)[1]
    end
end

# save data

h5open(data_folder_name*"results_.h5", "w") do file
    write(file, "results", res)
end

mkpath(output_folder_name* "/" * network * "/" *PN_name)
fig_file_name = output_folder_name * "grid_plot.pdf"
println("\nSaving grid plot in: $(fig_file_name)")

splot = heatmap(mu_vec, sigma_vec, res')
savefig(splot, fig_file_name)
