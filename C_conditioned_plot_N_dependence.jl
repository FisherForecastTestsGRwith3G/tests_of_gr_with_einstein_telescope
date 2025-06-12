using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# Fast modification of the original script C, to produce a plot of the the trend of the hierarchical upper limits as a function of the number of events, conditioned on the delta phi distribution.  

println("Running C_conditioned_plot_N_dependence.jl script: this script is a modified version of the C_conditioned_plot.jl script, not completely refined. You should look a the source code and set the correct parameters in the script (not just in the config files), to obtain the wanted result.")

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
min_number_of_realizations = 1
max_number_of_realizations = 40000
number_of_points_trend_plot = 20

# For simplicity, also to not clutter too much the plot, I will iterate over a single detector network
index_detector_network_to_use = 2

# Some settings, but OVERLOADED here in the code just below!!!
use_std_for_ribbon_in_log_space = configs["compute_mean_std_dev_in_log_space"]
use_std_for_ribbon_in_log_space = true # If true, the std will be computed in log space, otherwise in linear space. This is useful to avoid negative values in the ribbon (which lead to errors), but is slightly different from the standard deviation in linear space plotted in plot 6.
compute_mean_in_log_space = configs["compute_mean_std_dev_in_log_space"]
compute_mean_in_log_space = true # If true, the mean will be computed in log space, otherwise in linear space (as in plot 6, so default = false).
averageOverSeveralRealizations = configs["averageOverSeveralRealizations"]
averageOverSeveralRealizations = true
printEventsAsHorizontalLinesOrDensityPlot = false
n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = configs["n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events"]
n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = true # If true, the n_events and numberOfEventsSingleRealization refer directly to the observed events, and there is not poissonian/binomial noise due to drawing observed events from the catalog of all events. This is useful to remove the poissonian/binomial noise contribution to the shown error bars/ribbons.

#   Plot Settings (currently used for y labels only)
min_y = 10^-9
max_y = 10^2
set_y_lims = true
y_lims = (2e-8, 8e0)

# Compute the number of events used in the hierarchical analysis, conditioned on the delta phi distribution
# obtain a evenly spaced vector
n_events_realizations = Int.(floor.(exp.(log(min_number_of_realizations) .+ (log(max_number_of_realizations) - log(min_number_of_realizations)) .* (collect(1:number_of_points_trend_plot) .- 1) ./ (number_of_points_trend_plot - 1))))

# Remove duplicate values, if any
n_events_realizations = unique(n_events_realizations)
if length(n_events_realizations) < number_of_points_trend_plot
    println("Warning: automatically reduced the number of points in n_events_realizations (from $(number_of_points_trend_plot) to $(length(n_events_realizations))) by discarding duplicates values.")
    number_of_points_trend_plot = length(n_events_realizations)
end
println("Number of events used in the hierarchical analysis, conditioned on the delta phi distribution: $(n_events_realizations)")

# Plot some lines going as the square root of the number of events, to show the trend of the upper limits as a function of the number of events used in the hierarchical analysis, conditioned on the delta phi distribution.
# Here you can just indicate the different initial constants for these trend lines, such that they will be plotted accordingly
# just set initial_height_square_roots_trend_line = [] if you don't want to plot any trend line
# initial_height_square_roots_trend_line = (exp.(log(10^-5) .+ (log(10^-0) - log(10^-5)) .* (collect(1:6) .- 1) ./ (6 - 1))) .* sqrt.(n_events_realizations[1])
initial_height_square_roots_trend_line = (exp.(log(10^-9) .+ (log(10^4) - log(10^-9)) .* (collect(1:14) .- 1) ./ (14 - 1))) .* sqrt.(n_events_realizations[1])

# Create an empty array to store the results for the upper limits (mean and std for each of them... Eventually, if you have a single realization, the second parameter (std_dev) will be a NaN).
upperLimits = zeros(length(configs["network_list"]), number_of_points_trend_plot, length(configs["pn_waveforms"]), 4) # 4 = mean, std_dev, mean_in_log_space, std_dev_in_log_space
# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]

if configs["compute_mean_std_dev_in_log_space"]
    println("You are using the mean and std dev evaluated in log space!")
else
    println("You are using the mean and std dev evaluated in linear space!")
end

if !(n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events)
    println("Warning: this may be the only plot where you may want to set n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = true, in order to remove the poissonian/binomial noise contribution to the shown error bars/ribbons!")
end

# process all the plots
 for (index_nn, nn) in enumerate(configs["network_list"])

    println("\n"*"#"^81)
    println("Processing Network = $(nn)\n")

    if index_nn != index_detector_network_to_use
        println("Skipping network $(nn), as it is not the one selected for the trend plot.")
        continue
    end

    for index_realization in 1:number_of_points_trend_plot

        #Produce the plot evenly spaced in log scale        
        numberOfEventsSingleRealization = n_events_realizations[index_realization]

        println("Processing realization = $(index_realization) out of $(number_of_points_trend_plot) with $(numberOfEventsSingleRealization) events")

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
                n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events,
                printEventsAsHorizontalLinesOrDensityPlot = false, #I set this to false, as I will not plot the single events upper limits in this trend plot
                use_all_n_events_for_single_event_sample_distribution = false, # was configs["use_all_n_events_for_single_event_sample_distribution"], but this is irrelevant here
                print_info_catalog_realization = false # Skip printing the statistics, since this will be iterated many times, and you can just evaluate this as described below
            )

            # Skip printing the statistics, since this will be iterated many times, and you can just evaluate this as 
            # N_events_used = (N_events in catalog * probability_of_event_to_be_selected) +- (sqrt(N_events in catalog * probability_of_event_to_be_selected * (1 - probability_of_event_to_be_selected)))
            # since this is a binomial distribution (which for high N_events_used we may also approximate as a gaussian)
            # # Print statistics about the number_events_single_realization, if drawn from the catalog (since it induces a 'poissonian noise' - actually distributed as a binomial, given the high probability in ET - in the number of events per realization)
            # if !(configs["n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events"])
            #     println("Number of realizations: ", length(number_events_single_realization))
            #     println("Number of realization with non-zero number of events (and so used to obtain and plot the 90% upper bounds): ", sum(number_events_single_realization .> 0))
            #     println("Number of events requested per realization: ", numberOfEventsSingleRealizationUsed)
            #     println("Average number of events observed per realization: ", mean(number_events_single_realization))
            #     println("Standard deviation of the number of events per realization: ", std(number_events_single_realization))
            # end

            # Suppress the warning due to high numbre of iterations
            # if sum(number_events_single_realization .== 0) > 0
            #     println("Warning! There are $(sum(number_events_single_realization .== 0)) realizations with zero number of events, so they will not be used to obtain the upper bounds (so will not influence the average upper bound, nor its error bars)!")
            # end

            # I will only select the realization with a non zero number of events, so with a upper limit that is not a NaN
            vectorUpperLimits = vectorUpperLimits[.!isnan.(vectorUpperLimits)]

            # Evaluate the mean and std of the upper limits
            upperLimits[index_nn, index_realization, index_pno, :] = [mean(vectorUpperLimits), std(vectorUpperLimits), exp(mean(log.(vectorUpperLimits))), exp(std(log.(vectorUpperLimits)))]

            # Disable this
            # # Save the single events upper limits if required
            # if printEventsAsHorizontalLinesOrDensityPlot
            #     upperLimitSingleEvents[index_nn, index_pno, 1:length(upperLimitSingleEventsTemp)] = upperLimitSingleEventsTemp
            # end
        end
    end
end

# Call the specific plotting function

#Now I will plot upperLimits, for the single requested detector network, as a function of the number of events used in the hierarchical analysis, conditioned on the delta phi distribution.

# I setup a logscale-logscale plot, with the x-axis being the number of events used in the hierarchical analysis, and the y-axis being the upper limits.
set_common_plot_style()

# Settings:
plotHeight, plotWidth, plotDpi, padding = default_plot_dimensions() # Get the default plot dimensions

# Define a list of different markers and colors to enhance readability
markers, markersize, pointColors = get_markers_and_palette()


# Proper LaTeX labels
xlabel_str = L"N_{obs}"
ylabel_str = L"|\delta\varphi_{\!i}|" 
plotTitle = ""

# name_PN = collect(keys(pn_order_dic));
PN_orders = configs["pn_waveforms"];
name_PN = configs["pn_waveforms"];
network_names = configs["network_list"];
network_labels = labels_from_networks(network_names);
PN_labels = labels_from_PN_orders(PN_orders);
 
# Extend the arrays, eventually wrapping around
markers = markers[mod1.(1:length(PN_orders), length(markers))]
pointColors = pointColors[mod1.(1:length(PN_orders), length(pointColors))]

# ErrorBarsIntervalMultiplier increases the width of the error bars by this given amount 
# (the error bar here refers only to the estimation of the statistical scatter due to repeated different experiment realization, NOT to the probability threshold of the upper limit for the BeyondGR parameters, which you fixed in the C_conditioned_plot.jl code!)
# if ErrorBarsIntervalMultiplier = 1, then you are showing approx the 68% interval,
# if ErrorBarsIntervalMultiplier = 1.645, then you are showing approx the 90% interval
# if ErrorBarsIntervalMultiplier = 2, then you are showing approx the 95% interval

ErrorBarsIntervalMultiplier = 1.645;

final_plot = plot(
    xlabel=xlabel_str, 
    ylabel=ylabel_str,
    title=plotTitle,
    legend=:bottomright,
    #xticks=(1:length(pn_order_indices), PN_labels[pn_order_indices]),
    yscale=:log10, 
    xscale=:log10,
    xlims=(n_events_realizations[1], n_events_realizations[end]),  # Set x-axis limits to the first and last number of events used in the hierarchical analysis
    ylims=(set_y_lims ? y_lims : :auto),  # Set y-axis limits if requested   
    #left_margin = (index_subplot ==  1 ? small_margin_subplots.left_margin_with_ylabel : small_margin_subplots.left_margin),
    #right_margin = small_margin_subplots.right_margin,
    #top_margin = small_margin_subplots.top_margin,
    #bottom_margin = small_margin_subplots.bottom_margin,
    yticks=[10.0^i for i in floor(Int, log10(min_y)):ceil(Int, log10(max_y))],  # Set y-axis ticks for each 10^N value within the range
    size=(plotWidth, plotHeight), 
    padding = padding, 
    dpi=plotDpi,
    framestyle = :box,    
    minorgridcolor=:gray,  # Enable minor grid lines
    xminorgrid=true,  # Disable minor grid lines for the x-axis,     
    grid = true,
    gridalpha=0.4, 
    gridcolor=:gray,  # Set grid lines to be transparent or gray
    yminorgrid=true, 
    minorgridalpha=0.15
)

# Plot the trend lines for the square root of the number of events, if requested
if length(initial_height_square_roots_trend_line) > 0
    for (index_h, initial_height) in enumerate(initial_height_square_roots_trend_line)
        # Plot the square root trend line
        plot!(
            n_events_realizations, 
            initial_height ./ sqrt.(n_events_realizations), 
            #label="Trend line: " * string(initial_height) * " * sqrt(N)", 
            color=:gray, 
            alpha=0.5,
            label=(index_h == 1 ? L"\propto N_{obs}^{-1/2}" : ""),  # Only label the first trend line
            linestyle=:dash, 
            linewidth=2
        )
    end
end

# I plot the upper limits for each PN order, as a function of the number of events used in the hierarchical analysis, conditioned on the delta phi distribution.
for (index_pno, pno) in enumerate(PN_orders)
    # Extract the upper limits for the current PN order
    upperLimits_pno = upperLimits[index_detector_network_to_use, :, index_pno, (compute_mean_in_log_space ? 3 : 1)]  # Mean values
    upperLimits_std_pno = upperLimits[index_detector_network_to_use, :, index_pno, (use_std_for_ribbon_in_log_space ? 4 : 2)]  # Standard deviation values

    if use_std_for_ribbon_in_log_space
        lower = upperLimits_pno .- exp.(log.(upperLimits_pno) .- ErrorBarsIntervalMultiplier .* log.(upperLimits_std_pno))
        upper = exp.(log.(upperLimits_pno) .+ ErrorBarsIntervalMultiplier .* log.(upperLimits_std_pno)) .- upperLimits_pno
        ribbon_lower_upper_limits = hcat(lower, upper)
    else
        lower = upperLimits_pno .- ErrorBarsIntervalMultiplier .* upperLimits_std_pno
        upper = upperLimits_pno .+ ErrorBarsIntervalMultiplier .* upperLimits_std_pno
        ribbon_lower_upper_limits = hcat(lower, upper)
    end

    # Plot the mean upper limits with error bars (as a ribbon), but without markers
    plot!(
        n_events_realizations, 
        upperLimits_pno[:], 
        ribbon= (ribbon_lower_upper_limits[:, 1], ribbon_lower_upper_limits[:, 2]),
        label=PN_labels[index_pno], 
        color=pointColors[index_pno], 
        marker=:none,#markers[index_pno], 
        markersize=markersize, 
        fillalpha=0.2,
        linewidth=3, 
        linestyle=:solid
    )
end

# Save the combined plot to file
mkpath(output_folder_name)
println("Saving the conditioned upper limits plot to file: ", output_folder_name * "plot_trend_conditioned_delta_phi_upper_limits_" * simulation_tag * ".pdf")
savefig(final_plot, output_folder_name * "plot_trend_conditioned_delta_phi_upper_limits_" * simulation_tag * ".pdf")