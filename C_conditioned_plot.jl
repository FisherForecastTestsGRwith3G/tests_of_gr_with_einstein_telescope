using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# include required scripts 
#include("_hierarchical_dist.jl")
include("_hierarchical_deltaphi_dist.jl")
include("_parse_config.jl")
include("_conditioned_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForC(config_file_name)

# Fast modification of the script, to allow the possibility to plot several waveform models in the same plot (disguised as different detector networks)
overload_detector_networks_as_waveform_models = false # Default should be false! Set to true only if you know what you are doing!
if overload_detector_networks_as_waveform_models
    println("Warning: overloading the detector networks as waveform models! Usually this is not the standard behavior! Also, several config settings must be appopriately overloaded directly in the script! Set list_header_simulation_tags accordingly! Are you sure you want to proceed? (y/n)")
    # Ask for user input to confirm the change
    answer = readline()
    if answer != "y" && answer != "Y"
        println("Exiting the script!")
        exit(0)
    end

    # Settings to be correctly set in the script, if you want to overload the detector networks as waveform models
    # both waveform models must have data for the same detector networks! I assume so in the following code!
    list_header_simulation_tags = ["final_run_LVK", "final_run_LVK_PhenomD"] # As set in the config file, such that it indicates the correct folder where data is stored, list_header_simulation_tags = ["final_run_LVK", "final_run_LVK_PhenomD"]
    
    list_of_networks = [[x, y] for x in list_header_simulation_tags, y in configs["network_list"]]
    list_of_networks = vcat(list_of_networks...) 
else
    list_of_networks = configs["network_list"]
end

# Dictionaries to store the global indices for each network
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()

# Settings to perform the average over several realization of a given experiment, assuming a fixed number of observations for each run
averageOverSeveralRealizations = configs["averageOverSeveralRealizations"]
numberOfEventsSingleRealization = configs["numberOfEventsSingleRealization"] #Should be 12 ~ 15 to match LVK GWTC-3 Test of GR paper

println("Average over several realizations: $(averageOverSeveralRealizations)")
println("Number of events per realization: $(numberOfEventsSingleRealization)")

# Plot Settings
printEventsAsHorizontalLinesOrDensityPlot = true #I will do this only for the first realization, if there are more than 1 realizations
plotLVK_GWTC3_results = configs["LVK"] #Overlay the GWTC-3 results on the plot

# Create an empty array to store the results for the upper limits (mean and std for each of them... Eventually, if you have a single realization, the second parameter (std_dev) will be a NaN).
upperLimits = zeros(length(list_of_networks), length(configs["pn_waveforms"]), 4) # 4 = mean, std_dev, mean_in_log_space, std_dev_in_log_space
# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]
upperLimitSingleEvents = zeros(length(list_of_networks), length(configs["pn_waveforms"]), n_events_used)
#For the moment I will fill the array with NaN, so that if the number of events used is lower than n_events_used, I can simply discard the NaN
upperLimitSingleEvents .= NaN

if configs["compute_mean_std_dev_in_log_space"]
    println("You are using the mean and std dev evaluated in log space!")
else
    println("You are using the mean and std dev evaluated in linear space!")
end

# process all the plots
 for (index_nn, name_network) in enumerate(list_of_networks)

    println("\n"*"#"^81)
    println("Processing Network = $(name_network)\n")


    # set up folder names
    if overload_detector_networks_as_waveform_models # Overloading the detector networks as waveform models if requested
        data_folder_name = user_configs["path_output"]*name_network[1]*"/data/"
        output_folder_name =  user_configs["path_output"]*name_network[1]*"/plots/" 
        file_name = data_folder_name * name_network[2] *"/global_indices.h5"
        nn = name_network[2]
    else 
        data_folder_name = user_configs["path_output"]*simulation_tag*"/data/"
        output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"
        file_name = data_folder_name * name_network *"/global_indices.h5"
        nn = name_network
    end

    # get global index
    h5open(file_name, "r") do gi_file
        global_index_fisher[nn] = read(gi_file, "fisher")  
        global_index_snr[nn] = read(gi_file, "snr")
        global_index_total[nn] = read(gi_file, "total")
    end

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

        vectorUpperLimits, upperLimitSingleEventsTemp, number_events_single_realization, numberOfEventsSingleRealizationUsed = obtain_conditioned_upper_bounds(
            configs["n_events"], 
            global_index_total[nn], 
            data["dphi0_k"], 
            data["delta_k"], 
            averageOverSeveralRealizations = averageOverSeveralRealizations, 
            numberOfEventsSingleRealization = numberOfEventsSingleRealization, 
            n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = configs["n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events"],
            printEventsAsHorizontalLinesOrDensityPlot = printEventsAsHorizontalLinesOrDensityPlot,
            use_all_n_events_for_single_event_sample_distribution = configs["use_all_n_events_for_single_event_sample_distribution"]
        )

        # # Print statistics about the number_events_single_realization, if drawn from the catalog (since it induces a 'poissonian noise' - actually distributed as a binomial, given the high probability in ET - in the number of events per realization)            
        # # This should follow a binomial distribution (which for high N_events_used we may also approximate as a gaussian), with N_events_used = (N_events in catalog * probability_of_event_to_be_selected) +- (sqrt(N_events in catalog * probability_of_event_to_be_selected * (1 - probability_of_event_to_be_selected)))
        # if !(configs["n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events"])
        #     println("Number of realizations: ", length(number_events_single_realization))
        #     println("Number of realization with non-zero number of events (and so used to obtain and plot the 90% upper bounds): ", sum(number_events_single_realization .> 0))
        #     println("Number of events requested per realization: ", numberOfEventsSingleRealizationUsed)
        #     println("Average number of events observed per realization: ", mean(number_events_single_realization))
        #     println("Standard deviation of the number of events per realization: ", std(number_events_single_realization))
        # end

        if sum(number_events_single_realization .== 0) > 0
            println("Warning! There are $(sum(number_events_single_realization .== 0)) realizations with zero number of events, so they will not be used to obtain the upper bounds (so will not influence the average upper bound, nor its error bars)!")
        end

        # I will only select the realization with a non zero number of events, so with a upper limit that is not a NaN
        vectorUpperLimits = vectorUpperLimits[.!isnan.(vectorUpperLimits)]

        # Evaluate the mean and std of the upper limits
        upperLimits[index_nn, index_pno, :] = [mean(vectorUpperLimits), std(vectorUpperLimits), exp(mean(log.(vectorUpperLimits))), exp(std(log.(vectorUpperLimits)))]

        
        # Save the single events upper limits if required
        if printEventsAsHorizontalLinesOrDensityPlot
            upperLimitSingleEvents[index_nn, index_pno, 1:length(upperLimitSingleEventsTemp)] = upperLimitSingleEventsTemp
        end

    end
end

# Call the specific plotting function
# Pass the title as a LaTeXString, with L"\mathrm{Title\ text}"
title = configs["title"] == "" ? "" : latexstring(configs["title"])
final_conditioned_plot = plotConditionedUpperLimits(title, 
    upperLimits, 
    printEventsAsHorizontalLinesOrDensityPlot, 
    upperLimitSingleEvents, 
    plotLVK_GWTC3_results = plotLVK_GWTC3_results, 
    plot_samples_distribution = configs["plot_samples"], 
    compute_mean_std_dev_in_log_space = configs["compute_mean_std_dev_in_log_space"], 
    list_of_networks_names = (overload_detector_networks_as_waveform_models ? list_of_networks[:,1] .* "_" .* list_of_networks[:,2] : list_of_networks)
)

# Save the combined plot to file
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"
mkpath(output_folder_name)
println("Saving the conditioned upper limits plot to file: ", output_folder_name * "plot_conditioned_delta_phi_upper_limits_" * simulation_tag * ".pdf")
savefig(final_conditioned_plot, output_folder_name * "plot_conditioned_delta_phi_upper_limits_" * simulation_tag * ".pdf")