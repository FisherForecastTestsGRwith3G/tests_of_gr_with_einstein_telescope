using StatsPlots

# include required scripts 
include("_plot_style.jl")


#Results from figure 6 of https://arxiv.org/pdf/2112.06861
LVK_GWTC3_results = [0.75e-3, 0.06, 0.15, 0.1, 0.07, 0.55, 0.23, 0.48, 2.0, 1.0]

"""
The plotConditionedUpperLimits function plots the cumulative error on PN order for different networks and PN orders
plotTitle should be something like "Cumulative Error on PN Order - PhenomHM ET"
upperLimits should be a 3D array with dimensions (number of networks, number of PN orders, 2)
upperLimits[:, :, 1] should contain the mean values
upperLimits[:, :, 2] should contain the standard deviation values
printEventsAsHorizontalLinesOrDensityPlot should be a boolean indicating whether to plot individual events as horizontal lines or density plots
upperLimitSingleEvents should be a 3D array with dimensions (number of networks, number of PN orders, number of events)
plotLVK_GWTC3_results = True overlays the results for the LVK GWTC-3 results (yet it does not rescale the axis yet, so they may be out of the plotted region!)

In the config files,
"use_all_n_events_for_single_event_sample_distribution": allows to use (almost) all events for the single event sample distribution, instead of only the ones from a single realization [true/false]
"offset_x_axis_hierarchical_networks_upper_bounds": allows to offset the x-axis for the hierarchical upper bounds, so that they do not overlap if there are multiple networks [float, default 0.1]
"offset_x_axis_single_event_networks_upper_bounds": allows to offset the x-axis for the single events upper bounds (e.g. violin plots), so that they do not overlap if there are multiple networks [float, default 0.0],
"impose_y_axis_limits": allows to impose the y-axis limits for the plot, so that they are not automatically set based on the data [true/false]
"y_axis_limits": allows to set the y-axis limits for the plot, if impose_y_axis_limits is true [2x2 vector, outer vector represent the left and right subplots, the inner vector represents the lower and upper limits, e.g. [[1e-8, 1e-6],[2e-5, 1e-2]] ]
"violin_plots_bandwidth_std_multiplier": allows to set the bandwidth for the kernel density estimation used in the violin plots, as a multiplier of the standard deviation of the data [float, default 0.4]
"""
function plotConditionedUpperLimits(plotTitle, upperLimits, printEventsAsHorizontalLinesOrDensityPlot = false, upperLimitSingleEvents = nothing; plotLVK_GWTC3_results = false, plot_samples_distribution::Bool = false, plot_samples_min_n_events_for_violin_plots::Int64 = 20)

    # Settings:

    # ErrorBarsIntervalMultiplier increases the width of the error bars by this given amount 
    # (the error bar here refers only to the estimation of the statistical scatter due to repeated different experiment realization, NOT to the probability threshold of the upper limit for the BeyondGR parameters, which you fixed in the C_conditioned_plot.jl code!)
    # if ErrorBarsIntervalMultiplier = 1, then you are showing approx the 68% interval,
    # if ErrorBarsIntervalMultiplier = 1.645, then you are showing approx the 90% interval
    # if ErrorBarsIntervalMultiplier = 2, then you are showing approx the 95% interval
    
    set_common_plot_style()

    ErrorBarsIntervalMultiplier = 1.645;

    plotHeight, plotWidth, plotDpi, padding = default_plot_dimensions() # Get the default plot dimensions
    ratioFirstToTotalPlotMarginsIncluded = 1/4.;
    ratioFirstToTotalPlotPlotOnly = 0.125;
    thresholdNumberEventsAboveWhichToPlotDensity = plot_samples_min_n_events_for_violin_plots # Threshold for the number of events above which to plot density instead of horizontal lines
    plot_padding_top= 3;
    plot_padding_bottom = 4;
    alpha_level_single_events = 0.2;
    violin_plots_width = 0.5;  # Width of the violin plots

    # Define a list of different markers to enhance readability
    markers, markersize, pointColors = get_markers_and_palette()
    horizontalLineColorNetworks = pointColors

    # Proper LaTeX labels
    xlabel_str = "PN order"  # Use \mathrm for proper LaTeX rendering
    ylabel_str = L"|\delta\varphi_{\!i}|"


    println("\n"*"#"^81)
    println("Plotting the upper limits for the conditioned plot.")

    if printEventsAsHorizontalLinesOrDensityPlot && isnothing(upperLimitSingleEvents)
        throw(ArgumentError("upperLimitSingleEvents is nothing, but printEventsAsHorizontalLinesOrDensityPlot is true."))
        printEventsAsHorizontalLinesOrDensityPlot = false
    end

    # name_PN = collect(keys(pn_order_dic));
    PN_orders = configs["pn_waveforms"];
    network_names = configs["network_list"];
    network_labels = labels_from_networks(network_names);
    PN_labels = labels_from_PN_orders(PN_orders);

    # Extend the arrays, eventually wrapping around
    markers = markers[mod1.(1:length(network_names), length(markers))]
    pointColors = pointColors[mod1.(1:length(network_names), length(pointColors))]
    horizontalLineColorNetworks = horizontalLineColorNetworks[mod1.(1:length(network_names), length(horizontalLineColorNetworks))]

    # Calculate mean and standard deviation for each ii and jj
    # upperLimits = zeros(length(PN_orders), length(network_names))
    # std_cumulative_error = zeros(length(PN_orders), length(network_names))


    if size(upperLimits)[2] != 10
        throw(ArgumentError("You did not process all PN orders, for now this function cannot deal with results processed only partially!"))
    end

    #I check that there are no NaNs for the (eventually average) value of the upperLimits
    if any(isnan.(upperLimits[:, :, 1]))
        throw(ArgumentError("Warning: There are NaNs in the 90% upper limits (their value, eventually the average one)."))
    end


    # If all the standard deviations in upperLimits have been provided, i.e. there are no NaNs, then I will display the error bars as well!
    displayErrorBars = !any(isnan.(upperLimits[:, :, 2]))
    # If the standard deviation is NaN, then I will not display the error bars
    if displayErrorBars
        println("Displaying error bars")
    else
        println("Not displaying error bars, since some std dev values are NaN (probably you had only a single experiment realization, so you were not able to evaluate the std deviation)!")

        # I don't modify the code for the plots if displayErrorBars == false, since if the yerr is a NaN, the error bars are not displayed anyway!
        # Yet, if displayErrorBars == true, then I set all std dev to NaN, so that the error bars are not displayed
        upperLimits[:, :, 2] .= NaN
    end
    

    if configs["impose_y_axis_limits"]
        println("Imposing y-axis limits, using the provided values from the config file.")

        min_y = configs["y_axis_limits"][2][1]
        max_y = configs["y_axis_limits"][2][2]
        min_y_sub = configs["y_axis_limits"][1][1]
        max_y_sub = configs["y_axis_limits"][1][2]
    else
        println("Not imposing y-axis limits, using the data to set the limits automatically.")

        # Find the minimum and maximum values for y-axis limits, excluding the first point, and accounting for the padding due to the eventual error bars
        min_y = (minimum(upperLimits[:, 2:end, 1]) - ErrorBarsIntervalMultiplier * maximum(filter(!isnan, upperLimits[:, 2:end, 2]); init=0.0) ) / plot_padding_bottom  # Divide for some padding
        if min_y <= 0
            min_y = (minimum(upperLimits[:, 2:end, 1])) / plot_padding_bottom  # Divide for some padding
        end
        max_y = (maximum(upperLimits[:, 2:end, 1]) + ErrorBarsIntervalMultiplier * maximum(filter(!isnan, upperLimits[:, 2:end, 2]); init=0.0) ) * plot_padding_top  # Multiply for some padding

        # Find the minimum and maximum values for y-axis limits, excluding the first point
        min_y_sub = (minimum(upperLimits[:, 1,  1]) - ErrorBarsIntervalMultiplier * maximum(filter(!isnan, upperLimits[:, 1, 2]); init=0.0) ) / plot_padding_bottom  # Divide for some padding
        if min_y_sub <= 0
            min_y_sub = (minimum(upperLimits[:, 1, 1])) / plot_padding_bottom  # Divide for some padding
        end
        max_y_sub = (maximum(upperLimits[:, 1, 1]) + ErrorBarsIntervalMultiplier * maximum(filter(!isnan, upperLimits[:, 1, 2]); init=0.0) ) * plot_padding_top  # Multiply for some padding
    end

    # println("min_y = ", min_y)
    # println("max_y = ", max_y)
    # println("min_y_sub = ", min_y_sub)
    # println("max_y_sub = ", max_y_sub)
    # println("size(upperLimits) = ", size(upperLimits))


    
    # Initialize the main plot with white background and high resolution
    cc_main = plot(
        xlabel=xlabel_str, title=plotTitle,
        legend=:bottomright, xticks=(2:10, PN_labels[2:end]), 
        yscale=:log10, size=(plotWidth * (1. - ratioFirstToTotalPlotMarginsIncluded), plotHeight), dpi=plotDpi,
        xlims=(1.5 - maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]]), 10.5 + maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]])),
        grid=true, framestyle=:box,
        yticks=[10.0^i for i in floor(Int, log10(min_y)):ceil(Int, log10(max_y))],  # Set y-axis ticks for each 10^N value within the range
        ylims=(min_y, max_y),  # Automatically set the y-axis range
        gridalpha=0.5, gridcolor=:gray,  # Set grid lines to be transparent or gray
        yminorgrid=true, minorgridalpha=0.3, minorgridcolor=:gray,  # Enable minor grid lines
        xminorgrid=false  # Disable minor grid lines for the x-axis
    )

    # Initialize the subplot for the first point
    cc_sub = plot(
        xlabel=xlabel_str, 
        ylabel=ylabel_str,
        legend=false, xticks=([1], [PN_labels[1]]), 
        xlims=(0.5 - maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]]), 1.5 + maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]])),
        yscale=:log10, size=(plotWidth * ratioFirstToTotalPlotMarginsIncluded, plotHeight), dpi=plotDpi,
        grid=true, framestyle=:box,
        yticks=[10.0^i for i in floor(Int, log10(min_y_sub)) : ceil(Int, log10(max_y_sub))],  # Set y-axis ticks for the first point
        ylims=(min_y_sub, max_y_sub),  # Automatically set the y-axis range
        gridalpha=0.5, gridcolor=:gray,  # Set grid lines to be transparent or gray
        yminorgrid=true, minorgridalpha=0.3, minorgridcolor=:gray,  # Enable minor grid lines
        xminorgrid=false  # Disable minor grid lines for the x-axis
    )


    # Plot the single events as horizontal lines, if requested
    # Plotting them first, so they will be in the lowest plot layer


    if printEventsAsHorizontalLinesOrDensityPlot

        # Checks, if for any PN order and network, the number of events is above the Threshold, yet before I need to discard all NaNs!
        toPlotDensity = any(length(upperLimitSingleEvents[jj, ii, :][.!isnan.(upperLimitSingleEvents[jj, ii, :])]) > thresholdNumberEventsAboveWhichToPlotDensity for ii in 1:length(PN_orders) for jj in 1:length(network_names))

        for jj in 1:length(network_names)

            for ii in 1:1
                vecc = upperLimitSingleEvents[jj, ii, :]
                # I need to discard all NaNs!
                vecc = vecc[.!isnan.(vecc)]
                if plot_samples_distribution == true
                    if !toPlotDensity
                        # Plot the single events as horizontal lines. If required, I also slightly offset the horizontal lines to avoid overlapping with the main plot

                        scatter!(
                            cc_sub, fill(( isnothing(configs["offset_x_axis_single_event_networks_upper_bounds"]) ? ii : ii + configs["offset_x_axis_single_event_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))) , length(vecc)), vecc, label="",
                            marker=:hline, markersize=5, markerstrokewidth=1, alpha=alpha_level_single_events, color=horizontalLineColorNetworks[jj]
                        )
                    else
                        # Plot violin plots to show the density distribution of the events, as a function of the upper limits (along the vertical axis), instead of plotting each single event
                        # Actually, since the violin plots do not work well on a log scale, I plot them on a separate linear y axis, but where I perform myself the transformation to the log space
                        # Add a second y-axis (linear, no label, no grid, no ticks)
                        cc_sub2 = twinx(cc_sub)
                        xlims!(cc_sub2, (0.5 - maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]]), 1.5 + maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]])))
                        ylims!(cc_sub2, log.((min_y_sub, max_y_sub)))
                        Plots.xlabel!(cc_sub2, "")

                        violin!(
                            cc_sub2, 
                            fill((ii + configs["offset_x_axis_single_event_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))), length(vecc)), 
                            log.(vecc), 
                            ylabel = "",
                            yticks = [],
                            #grid = false,
                            #legend = false,
                            #framestyle = :none,
                            orientation=:vertical, 
                            width=violin_plots_width, 
                            alpha=alpha_level_single_events, 
                            color=horizontalLineColorNetworks[jj], 
                            label="", 
                            linewidth = 0.0, 
                            bandwidth=std(log.(vecc)) * configs["violin_plots_bandwidth_std_multiplier"]
                        )
                    end
                end
            end
        
            for ii in 2:length(PN_orders)
                vecc = upperLimitSingleEvents[jj, ii, :]
                # I need to discard all NaNs!
                vecc = vecc[.!isnan.(vecc)]
                if plot_samples_distribution == true

                    if !toPlotDensity
                        # Plot the single events as horizontal lines
                        scatter!(
                            cc_main, fill(( isnothing(configs["offset_x_axis_single_event_networks_upper_bounds"]) ? ii : ii + configs["offset_x_axis_single_event_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))), length(vecc)), vecc, label="",
                            marker=:hline, markersize=5, markerstrokewidth=1, alpha=alpha_level_single_events, color=horizontalLineColorNetworks[jj]
                        )
                    else
                        # Plot the density of the single events
                        # Actually, since the violin plots do not work well on a log scale, I plot them on a separate linear y axis, but where I perform myself the transformation to the log space
                        # Add a second y-axis (linear, no label, no grid, no ticks)
                        cc_main2 = twinx(cc_main)
                        xlims!(cc_main2, (1.5 - maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]]), 10.5 + maximum([configs["offset_x_axis_single_event_networks_upper_bounds"], configs["offset_x_axis_hierarchical_networks_upper_bounds"]])))
                        ylims!(cc_main2, log.((min_y, max_y)))
                        Plots.xlabel!(cc_main2, "")

                        violin!(
                            cc_main2, 
                            fill((ii + configs["offset_x_axis_single_event_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))), length(vecc)), 
                            log.(vecc),
                            ylabel = "",
                            yticks = [],
                            #grid = false,
                            #legend = false,
                            #framestyle = :none,
                            orientation=:vertical, 
                            width=violin_plots_width, 
                            alpha=alpha_level_single_events, 
                            color=horizontalLineColorNetworks[jj], 
                            label="", 
                            linewidth = 0.0, 
                            bandwidth=std(log.(vecc)) * configs["violin_plots_bandwidth_std_multiplier"]
                        )
                    end
                end
            end

            #if !toPlotDensity
            if plot_samples_distribution == true

                scatter!(
                    cc_main, [NaN], [NaN],
                    marker=:hline, markersize=5, markerstrokewidth=(toPlotDensity ? 20 : 1), alpha=alpha_level_single_events, color=horizontalLineColorNetworks[jj], label="Single event bounds ("*network_labels[jj]*")"
                )
            end
            #else
            #   violin!(cc_main, [NaN, NaN], [NaN, NaN], orientation=:vertical, width=violin_plots_width, alpha=alpha_level_single_events, color=horizontalLineColorNetworks[jj], label="Single event bounds ("*network_names[jj]*")" )
            #end
            
        end

        
        
    end


    # Plot the first point in the subplot
    for jj in 1:length(network_names)
        scatter!(
            cc_sub, [1] .+ ( isnothing(configs["offset_x_axis_hierarchical_networks_upper_bounds"]) ? 0 : configs["offset_x_axis_hierarchical_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))), [upperLimits[jj, 1, 1]],
            yerr= ErrorBarsIntervalMultiplier * [upperLimits[jj, 1, 2]],  # Add error bars
            label=network_labels[jj],
            marker=markers[jj],  # Cycle through marker list
            color=pointColors[jj],
            markersize=markersize, markerstrokewidth=1, markerstrokecolor=:black,
            alpha=0.9  # Slight transparency for overlapping points
        )
    end

    # Plot the remaining points in the main plot
    for jj in 1:length(network_names)
        scatter!(
            cc_main, (2:length(PN_orders)) .+ ( isnothing(configs["offset_x_axis_hierarchical_networks_upper_bounds"]) ? 0 : configs["offset_x_axis_hierarchical_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names))), upperLimits[jj, 2:end, 1],
            yerr= ErrorBarsIntervalMultiplier * upperLimits[jj, 2:end, 2],  # Add error bars
            label=network_labels[jj],
            marker=markers[jj],  # Cycle through marker list
            color=pointColors[jj],
            markersize=markersize, markerstrokewidth=1, markerstrokecolor=:black,
            alpha=0.9  # Slight transparency for overlapping points
        )
    end



    
    if(plotLVK_GWTC3_results)
        # Plot the LVK GWTC-3 Results
        # Add a dummy plot for the legend entry in the main plot

        # Plot the GWTC-3 results in the main plot
        scatter!(
            cc_main, 2:length(PN_orders), LVK_GWTC3_results[2:end],
            label="LVK GWTC-3", marker=:diamond, markersize=markersize, markerstrokewidth=1, markerstrokecolor=:black, color=:darkblue
        )

        # Plot the first GWTC-3 result in the subplot
        scatter!(
            cc_sub, [1], [LVK_GWTC3_results[1]],
            label="LVK GWTC-3", marker=:diamond, markersize=markersize, markerstrokewidth=1, markerstrokecolor=:black, color=:darkblue
        )
    end


    # Combine the main plot and the subplot
    final_plot = plot(cc_sub, cc_main, layout = @layout([grid(1, 2, widths = [ratioFirstToTotalPlotPlotOnly, 1 - ratioFirstToTotalPlotPlotOnly])])) #, top_margin=2mm, bottom_margin=2mm, left_margin=2mm, right_margin=2mm)
    plot!(
        final_plot, 
        title=plotTitle, 
        #xlabel=xlabel_str, 
        size=(plotWidth, plotHeight), 
        padding = padding, 
        dpi=plotDpi,
        framestyle = :box,
        grid = true,
        gridalpha=0.5, 
        gridcolor=:gray,  # Set grid lines to be transparent or gray
        yminorgrid=true, 
        minorgridalpha=0.3
    ) #, ylabel=ylabel_str, top_margin=2mm, bottom_margin=2mm, left_margin=2mm, right_margin=2mm)

    # Return the final plot
    return final_plot
end
