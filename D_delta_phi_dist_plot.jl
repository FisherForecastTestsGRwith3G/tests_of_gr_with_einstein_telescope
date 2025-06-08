using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing
using Base.Threads

# include required scripts 
include("_parse_config.jl")
include("_hierarchical_deltaphi_dist.jl")
include("_delta_phi_dist_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 1
    config_file_name = ARGS[2]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForD(config_file_name)

# check if we should run the MCMC posterior sampling or load the results from disk
perform_MCMC_sampling = true
if length(ARGS) > 0
    if ARGS[1] == "1"
        perform_MCMC_sampling = true
    elseif ARGS[1] == "0"
        perform_MCMC_sampling = false
    else
        @warn "There was a first input argument handed. But it was $(ARGS[1]) and not \"1\" or \"0\" and thus ignored"
    end
end

if perform_MCMC_sampling
        println("MCMC sampling of the posterior will be performed!")
else
        println("MCMC sampling of the posterior will not be performed, and instead will be loaded from disk (saves time if you already evaluated them before)!")
end

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

# process all the networks and PN orders
for (index_nn, nn) in enumerate(configs["network_list"])

    # Filename to save or load the posterior distribution for delta_phi to/from disk
    filename_samples = data_folder_name * "script_D/" * nn * "/posterior_distributions_deltaphi.h5"

    if perform_MCMC_sampling
        println("\n"*"#"^81)
        println("Processing Network = $(nn)\n")
        
        @time begin
            @threads for index_pno = 1:length(configs["pn_waveforms"])
                pno = configs["pn_waveforms"][index_pno]
                
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

        println("Saving the posterior distribution for delta_phi to disk, in $(filename_samples)")
        mkpath(data_folder_name * "script_D/" * nn  * "/")
        # save to .h5 file 
        h5open(filename_samples, "w") do file
            write(file, "samples_posterior_dist_deltaphi", posterior_dist_deltaphi)
            if configs["plot_conditioned_distribution"]
                write(file, "samples_posterior_dist_deltaphi_conditioned", posterior_dist_deltaphi_conditioned)
            end
        end
    else
        # If we do not perform the MCMC sampling, I will load the results from disk
        println("Loading the posterior distribution for delta_phi from disk, from file $(filename_samples)")
        global posterior_dist_deltaphi, posterior_dist_deltaphi_conditioned = h5open(filename_samples, "r") do file
            posterior_dist_deltaphi = read(file, "samples_posterior_dist_deltaphi")
            if configs["plot_conditioned_distribution"]
                posterior_dist_deltaphi_conditioned = read(file, "samples_posterior_dist_deltaphi_conditioned")
            else
                posterior_dist_deltaphi_conditioned = nothing
            end
            return posterior_dist_deltaphi, posterior_dist_deltaphi_conditioned
        end
        # Check if the loaded posterior distribution has the correct dimensions
        if size(posterior_dist_deltaphi, 1) != length(configs["network_list"]) || size(posterior_dist_deltaphi, 2) != length(configs["pn_waveforms"]) || size(posterior_dist_deltaphi, 3) != configs["mcmc_samples"]
            throw(ArgumentError("The loaded posterior distribution for delta_phi has incorrect dimensions. Expected $(length(configs["network_list"])), $(length(configs["pn_waveforms"])), $(configs["mcmc_samples"]), but got $(size(posterior_dist_deltaphi))."))
        end
        if configs["plot_conditioned_distribution"] && (posterior_dist_deltaphi_conditioned === nothing || size(posterior_dist_deltaphi_conditioned, 1) != length(configs["network_list"]) || size(posterior_dist_deltaphi_conditioned, 2) != length(configs["pn_waveforms"]) || size(posterior_dist_deltaphi_conditioned, 3) != configs["mcmc_samples"])
            throw(ArgumentError("The loaded posterior distribution for delta_phi conditioned has incorrect dimensions. Expected $(length(configs["network_list"])), $(length(configs["pn_waveforms"])), $(configs["mcmc_samples"]), but got $(size(posterior_dist_deltaphi_conditioned))."))
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