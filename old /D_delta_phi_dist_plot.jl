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
number_of_expected_realizations = (configs["evaluate_different_realizations"] ? Int(floor(configs["n_events"] / configs["number_of_events_single_realization"])) : 1)
# Samples from the posterior distribution of delta_phi
posterior_dist_deltaphi = zeros(number_of_expected_realizations, length(configs["network_list"]), length(configs["pn_waveforms"]), configs["mcmc_samples"])

if configs["plot_conditioned_distribution"]
    # Samples from the posterior distribution of delta_phi, conditioned on sigma = 0 for the hyperparameter distribution. 
    # In this case I know the analytical form of the pdf, but for simplicity I will just sample from it to produce the violin plots.
    posterior_dist_deltaphi_conditioned = zeros(number_of_expected_realizations, length(configs["network_list"]), length(configs["pn_waveforms"]), configs["mcmc_samples"])
else
    posterior_dist_deltaphi_conditioned = nothing
end

# process all the networks and PN orders
for (index_nn, nn) in enumerate(configs["network_list"])

    # Filename to save or load the posterior distribution for delta_phi to/from disk

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

                dphi0_k_realization, delta_k_realization, numberOfRealization, number_events_single_realization = obtain_dphi0k_deltak_from_realizations(
                    configs["n_events"], 
                    global_index_total[nn], 
                    data["dphi0_k"], 
                    data["delta_k"], 
                    averageOverSeveralRealizations = configs["evaluate_different_realizations"], 
                    numberOfEventsSingleRealization = configs["number_of_events_single_realization"],
                    n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = configs["n_events_and_number_of_events_single_realization_refer_directly_to_observed_events"],
                    print_info_catalog_realization = true
                )

                if numberOfRealization > number_of_expected_realizations
                    # If we have more realizations than expected, we will only use the first number_of_expected_realizations realizations
                    println("Warning: There are $(numberOfRealization) realizations available, but only $(number_of_expected_realizations) were expected (and array space was allocated for). Using only the first $(number_of_expected_realizations) realizations to avoid out-of-bounds errors.")
                    numberOfRealization = number_of_expected_realizations
                end

                for realization_index in 1:numberOfRealization
                    # Now I iterate over each realization

                    if number_events_single_realization[realization_index] == 0
                        # No events in realization, so I will skip this realization and set the upper limit to NaN
                        println("Realization $(realization_index) has no events, this may be problematic! Skipping this realization, while setting its entries to NaNs!.")
                        posterior_dist_deltaphi[realization_index, index_nn, index_pno, :] .= NaN
                        if configs["plot_conditioned_distribution"]
                            posterior_dist_deltaphi_conditioned[realization_index, index_nn, index_pno, :] .= NaN
                        end
                        continue
                    end

                    # Obtain samples for posterior distributions for delta_phi
                    posterior_dist_deltaphi[realization_index, index_nn, index_pno, :] = obtain_samples_delta_phi_pdf(dphi0_k_realization[realization_index], delta_k_realization[realization_index]; n_samples = configs["mcmc_samples"], burn_in = configs["mcmc_burnin"], debug_folder_name = debug_folder_name * nn * "/realization_" * string(realization_index) * "/", pnorder = pn_order_dic[pno][2])

                    if configs["plot_conditioned_distribution"]
                        # Obtain samples for posterior distributions for delta_phi, condition on sigma = 0 in the hierarchical distribution
                        posterior_dist_deltaphi_conditioned[realization_index, index_nn, index_pno, :] = obtain_samples_delta_phi_pdf_conditioned(dphi0_k_realization[realization_index], delta_k_realization[realization_index]; n_samples = configs["mcmc_samples"])
                    end

                end



            end
        end

        for realization_index in 1:number_of_expected_realizations
            filename_samples = data_folder_name * "script_D/" * nn  * "/realization_" * string(realization_index) * "/posterior_distributions_deltaphi.h5"
            println("Saving the posterior distribution for delta_phi to disk, in $(filename_samples)")
            mkpath(data_folder_name * "script_D/" * nn  * "/realization_" * string(realization_index) * "/")
            # save to .h5 file 
            h5open(filename_samples, "w") do file
                write(file, "samples_posterior_dist_deltaphi", posterior_dist_deltaphi[realization_index, index_nn, :, :])
                if configs["plot_conditioned_distribution"]
                    write(file, "samples_posterior_dist_deltaphi_conditioned", posterior_dist_deltaphi_conditioned[realization_index, index_nn, :, :])
                end
            end
        end

    else
        for realization_index in 1:number_of_expected_realizations
            filename_samples = data_folder_name * "script_D/" * nn  * "/realization_" * string(realization_index) * "/posterior_distributions_deltaphi.h5"
            # If we do not perform the MCMC sampling, I will load the results from disk
            println("Loading the posterior distribution for delta_phi from disk, from file $(filename_samples)")
            # Check if the file exists
            if !isfile(filename_samples)
                throw(ArgumentError("The file $(filename_samples) does not exist. Please run the MCMC sampling first, or check the file path, or check that the current expected number of realizations ($(number_of_expected_realizations)) matches the precomputed ones."))
            end
            local network_posterior_dist_deltaphi, network_posterior_dist_deltaphi_conditioned = h5open(filename_samples, "r") do file
                network_posterior_dist_deltaphi = read(file, "samples_posterior_dist_deltaphi")
                if configs["plot_conditioned_distribution"]
                    network_posterior_dist_deltaphi_conditioned = read(file, "samples_posterior_dist_deltaphi_conditioned")
                else
                    network_posterior_dist_deltaphi_conditioned = nothing
                end
                return network_posterior_dist_deltaphi, network_posterior_dist_deltaphi_conditioned
            end
            # Check if the loaded posterior distribution has the correct dimensions
            if size(network_posterior_dist_deltaphi, 1) != length(configs["pn_waveforms"]) || size(network_posterior_dist_deltaphi, 2) != configs["mcmc_samples"]
                throw(ArgumentError("The loaded posterior distribution for delta_phi has incorrect dimensions. Expected $(length(configs["pn_waveforms"])), $(configs["mcmc_samples"]), but got $(size(network_posterior_dist_deltaphi))."))
            end
            if configs["plot_conditioned_distribution"] && (posterior_dist_deltaphi_conditioned === nothing || size(network_posterior_dist_deltaphi_conditioned, 1) != length(configs["pn_waveforms"]) || size(network_posterior_dist_deltaphi_conditioned, 2) != configs["mcmc_samples"])
                throw(ArgumentError("The loaded posterior distribution for delta_phi conditioned has incorrect dimensions. Expected $(length(configs["pn_waveforms"])), $(configs["mcmc_samples"]), but got $(size(network_posterior_dist_deltaphi_conditioned))."))
            end
            # Assign the loaded posterior distributions to the main variable
            posterior_dist_deltaphi[realization_index, index_nn, :, :] = network_posterior_dist_deltaphi
            if configs["plot_conditioned_distribution"]
                posterior_dist_deltaphi_conditioned[realization_index, index_nn, :, :] = network_posterior_dist_deltaphi_conditioned
            end
        end
    end
end

mkpath(output_folder_name)
# Call the specific plotting function
if configs["average_posteriors_from_different_realizations"]
    # Pass the title as a LaTeXString, with L"\mathrm{Title\ text}"
    local title = configs["title"] == "" ? "" : latexstring(configs["title"])

    # Here, to perform the average, I translate each MCMC chain by its mean value, and then sum the different chains... most surely this is not the best way to do it, but it is the simplest one 
    # (yet we lose information about the variance of the position of the mean, even thought it should be also somehow encoded in the spread of the posterior itself, since we expect our results to not be biased)
    # TODO: find a better way to average the posterior distributions

    means_posterior_dist_deltaphi = mean(posterior_dist_deltaphi, dims = 4) # average over the MCMC samples
    translated_posterior_dist_deltaphi = posterior_dist_deltaphi .- means_posterior_dist_deltaphi
    # Concatenate the first and fourth dimensions into a new third dimension. Now, averaged_posterior_dist_deltaphi has size (N2, N3, N1*N4)
    averaged_posterior_dist_deltaphi = reshape(permutedims(translated_posterior_dist_deltaphi, (2, 3, 1, 4)), size(translated_posterior_dist_deltaphi, 2), size(translated_posterior_dist_deltaphi, 3), :)
    # Print the average value of the means, and their std dev
    println("Average value of the means of the posterior distributions for delta_phi (over networks, PN orders): ", dropdims(mean(means_posterior_dist_deltaphi, dims=(1,4)),dims=(1,4)))
    println("Standard deviation of the means of the posterior distributions for delta_phi (over networks, PN orders): ", dropdims(std(means_posterior_dist_deltaphi, dims=(1,4)),dims=(1,4)))

    if configs["plot_conditioned_distribution"]
        # For the conditioned posterior distribution, I will do the same
        means_posterior_dist_deltaphi_conditioned = mean(posterior_dist_deltaphi_conditioned, dims = 4) # average over the MCMC samples
        translated_posterior_dist_deltaphi_conditioned = posterior_dist_deltaphi_conditioned .- means_posterior_dist_deltaphi_conditioned
        averaged_posterior_dist_deltaphi_conditioned = reshape(permutedims(translated_posterior_dist_deltaphi_conditioned, (2, 3, 1, 4)), size(translated_posterior_dist_deltaphi_conditioned, 2), size(translated_posterior_dist_deltaphi_conditioned, 3), :)
        println("Average value of the means of the posterior distributions for delta_phi conditioned (over networks, PN orders): ", dropdims(mean(means_posterior_dist_deltaphi_conditioned, dims=(1,4)),dims=(1,4)))
        println("Standard deviation of the means of the posterior distributions for delta_phi conditioned (over networks, PN orders): ", dropdims(std(means_posterior_dist_deltaphi_conditioned, dims=(1,4)),dims=(1,4)))
    end

    local delta_phi_posterior_plot = plotDeltaPhiPosterior(title, averaged_posterior_dist_deltaphi, plot_conditioned_distribution = configs["plot_conditioned_distribution"], posterior_dist_deltaphi_conditioned = averaged_posterior_dist_deltaphi_conditioned, subplots_pn_order_grouping = (configs["subplots_pn_order_grouping"] === nothing ? nothing : [Vector{Int}(x) for x in configs["subplots_pn_order_grouping"]]))

    # Save the combined plot to file
    println("Saving the plot to: ", output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * ".pdf")
    savefig(delta_phi_posterior_plot, output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * ".pdf")
else
    # I will plot all the posterior distributions separately, for each realization
    for realization_index in 1:number_of_expected_realizations
        # Pass the title as a LaTeXString, with L"\mathrm{Title\ text}"
        local title = configs["title"] == "" ? "" : latexstring(configs["title"])
        local delta_phi_posterior_plot = plotDeltaPhiPosterior(title, posterior_dist_deltaphi[realization_index, :, :, :], plot_conditioned_distribution = configs["plot_conditioned_distribution"], posterior_dist_deltaphi_conditioned = (isnothing(posterior_dist_deltaphi_conditioned) ? posterior_dist_deltaphi_conditioned : posterior_dist_deltaphi_conditioned[realization_index, :, :, :]), subplots_pn_order_grouping = (configs["subplots_pn_order_grouping"] === nothing ? nothing : [Vector{Int}(x) for x in configs["subplots_pn_order_grouping"]]))

        # Save the combined plot to file
        println("Saving the plot to: ", output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * "_realization_" * string(realization_index) * ".pdf")
        savefig(delta_phi_posterior_plot, output_folder_name * "plot_delta_phi_posterior_dist_" * simulation_tag * "_realization_" * string(realization_index) * ".pdf")
    end
end
