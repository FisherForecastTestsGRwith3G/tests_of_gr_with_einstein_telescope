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

# set up folder names
data_folder_name = user_configs["path_output"]*simulation_tag*"/data/"
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"

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


# Settings to perform the average over several realization of a given experiment, assuming a fixed number of observations for each run
averageOverSeveralRealizations = configs["averageOverSeveralRealizations"]
numberOfEventsSingleRealization = configs["numberOfEventsSingleRealization"] #Should be 12 ~ 15 to match LVK GWTC-3 Test of GR paper

println("Average over several realizations: $(averageOverSeveralRealizations)")
println("Number of events per realization: $(numberOfEventsSingleRealization)")

# Plot Settings
printEventsAsHorizontalLinesOrDensityPlot = true #I will do this only for the first realization, if there are more than 1 realizations
plotLVK_GWTC3_results = configs["LVK"] #Overlay the GWTC-3 results on the plot

# Create an empty array to store the results for the upper limits (mean and std for each of them... Eventually, if you have a single realization, the second parameter (std_dev) will be a NaN).
upperLimits = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), 2)
# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]
upperLimitSingleEvents = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), n_events_used)
#For the moment I will fill the array with NaN, so that if the number of events used is lower than n_events_used, I can simply discard the NaN
upperLimitSingleEvents .= NaN

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
        if averageOverSeveralRealizations
            numberOfRealization = Int(floor(n_events_used / numberOfEventsSingleRealization))
            n_events_used = numberOfEventsSingleRealization * numberOfRealization
            dphi0_k = reshape(dphi0_k[1:n_events_used], numberOfRealization,  numberOfEventsSingleRealization)
            delta_k = reshape(delta_k[1:n_events_used], numberOfRealization, numberOfEventsSingleRealization)

            println("Splitting the dataset into $(numberOfRealization) realizations of $(numberOfEventsSingleRealization) events")
            println("Number of events used effectively: $(n_events_used)")
            
        else
            numberOfRealization = 1
            dphi0_k = reshape(dphi0_k[1:n_events_used], numberOfRealization,  n_events_used)
            delta_k = reshape(delta_k[1:n_events_used], numberOfRealization, n_events_used)
        end

        vectorUpperLimits = zeros(numberOfRealization)
        for realization_index in 1:numberOfRealization
             # Access the realization for dphi0_k_realization
            dphi0_k_realization = dphi0_k[realization_index, :]
            
            # Access the realization for delta_k_realization
            delta_k_realization = delta_k[realization_index, :]

            #Evaluate the 90% upper limit given this single realization
            vectorUpperLimits[realization_index] = get90PctUpperLimitMuDistTIGER(dphi0_k_realization, delta_k_realization, 0.9)

            #Save the single events 90% upper limits for a given realization, if requested
            if printEventsAsHorizontalLinesOrDensityPlot && realization_index == 1
                # The hierarchical analysis, conditioned on sigma = 0, implemented in getMuDistTIGER, should work just fine even if working with a single event
                upperLimitSingleEventsTemp = map((x, y) -> get90PctUpperLimitMuDistTIGER([x], [y], 0.9), dphi0_k_realization, delta_k_realization)
                upperLimitSingleEvents[index_nn, index_pno, 1:length(upperLimitSingleEventsTemp)] = upperLimitSingleEventsTemp
            end

        end
        # Evaluate the mean and std of the upper limits
        upperLimits[index_nn, index_pno, :] = [mean(vectorUpperLimits), std(vectorUpperLimits)]

    end
end

# Call the specific plotting function
# Pass the title as a LaTeXString, with L"\mathrm{Title\ text}"
title = configs["title"] == "" ? "" : latexstring(configs["title"])
final_conditioned_plot = plotConditionedUpperLimits(title, upperLimits, printEventsAsHorizontalLinesOrDensityPlot, upperLimitSingleEvents, plotLVK_GWTC3_results = plotLVK_GWTC3_results)

# Save the combined plot to file
mkpath(output_folder_name)
savefig(final_conditioned_plot, output_folder_name * "plot_conditioned_upper_limits_" * simulation_tag * ".pdf")