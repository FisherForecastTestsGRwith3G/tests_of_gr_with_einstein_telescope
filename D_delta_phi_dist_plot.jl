using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# include required scripts 
include("_parse_config.jl")
include("_hierarchical_deltaphi_dist.jl")
include("_delta_phi_dist_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForD(config_file_name)

# set up folder names
data_folder_name = user_configs["path_output"]*simulation_tag*"/data/"
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"
debug_folder_name =  user_configs["path_output"]*simulation_tag*"/debug/script_D/"

# get global index
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()
for nn in configs["network_list"]
    file_name = data_folder_name * nn *"/global_indices.h5"
    h5open(file_name, "r") do gi_file
        global_index_fisher[nn] = read(gi_file, "fisher")  
        global_index_snr[nn] = read(gi_file, "snr")
        global_index_total[nn] = read(gi_file, "total")
    end
end

# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]
# Samples from the posterior distribution of delta_phi
posterior_dist_deltaphi = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), configs["mcmc_samples"])

if configs["plot_conditioned_distribution"]
    # Samples from the posterior distribution of delta_phi, conditioned on sigma = 0 for the hyperparameter distribution. 
    # In this case I know the analytical form of the pdf, but for simplicity I will just sample from it to produce the violin plots.
    posterior_dist_deltaphi_conditioned = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), configs["mcmc_samples"])
else
    posterior_dist_deltaphi_conditioned = nothing
end

# process all the plots
 for (index_nn, nn) in enumerate(configs["network_list"])

    println("\n"*"#"^81)
    println("Processing Network = $(nn)\n")

    for (index_pno, pno) in enumerate(configs["pn_waveforms"])

        println("\n"*"#"^81)
        println("Processing PN = $(pno)\n")
        
        # load the data
        folder = data_folder_name * nn * "/pn_" * pn_order_dic[pno][2]
        filename = folder * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # extract only valid data-points via global index

        local n_events_used = configs["n_events"]
        if n_events_used > sum(global_index_total[nn])
            println("Warning: There are at most $(sum(global_index_total[nn]))-datapoints available")
            n_events_used = sum(global_index_total[nn])
        else
            println("Events used: $(n_events_used)")
        end

        dphi0_k = data["dphi0_k"][global_index_total[nn]]
        dphi0_k = dphi0_k[1:n_events_used] 
        delta_k = data["delta_k"][global_index_total[nn]]
        delta_k = delta_k[1:n_events_used]

        # Obtain samples for posterior distributions for delta_phi
        posterior_dist_deltaphi[index_nn, index_pno, :] = obtain_samples_delta_phi_pdf(dphi0_k, delta_k; n_samples = configs["mcmc_samples"], burn_in = configs["mcmc_burnin"], debug_folder_name = debug_folder_name * nn * "/", pnorder = pn_order_dic[pno][2])

        if configs["plot_conditioned_distribution"]
            # Obtain samples for posterior distributions for delta_phi, condition on sigma = 0 in the hierarchical distribution
            posterior_dist_deltaphi_conditioned[index_nn, index_pno, :] = obtain_samples_delta_phi_pdf_conditioned(dphi0_k, delta_k; n_samples = configs["mcmc_samples"])
        end

    end
end

# Call the specific plotting function
# Pass the title as a LaTeXString, with L"\mathrm{Title\ text}"
title = configs["title"] == "" ? "" : latexstring(configs["title"])
delta_phi_posterior_plot = plotDeltaPhiPosterior(title, posterior_dist_deltaphi, plot_conditioned_distribution = configs["plot_conditioned_distribution"], posterior_dist_deltaphi_conditioned = posterior_dist_deltaphi_conditioned, subplots_pn_order_grouping = (configs["subplots_pn_order_grouping"] === nothing ? nothing : [Vector{Int}(x) for x in configs["subplots_pn_order_grouping"]]))

# Save the combined plot to file
mkpath(output_folder_name)
println("Saving the plot to: ", output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * ".pdf")
savefig(delta_phi_posterior_plot, output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * ".pdf")