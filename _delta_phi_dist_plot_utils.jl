using StatsPlots

# include required scripts 
include("_plot_style.jl")

"""
The plotDeltaPhiPosterior function plots the posterior distributions ('marginalized' and eventually 'conditioned') for the PN deformation coefficients delta_phi, for several networks and PN orders
- plotTitle should be something like "Posterior distributions for delta varphi"
- posterior_dist_deltaphi should be a 3D array with dimensions (number of networks, number of PN orders, number of samples)
- plot_conditioned_distribution is a float that indicates whether to plot the conditioned distribution (default is false)
- posterior_dist_deltaphi_conditioned should be either nothing, if plot_conditioned_distribution == false, or a 3D array with dimensions (number of networks, number of PN orders, number of samples) otherwise
- subplots_pn_order_grouping should be either nothing (for automatically plotting all provided PN orders in the same subplot), or a 2D vector if you want to plot all different PN orders in different subplots (in that case, each subvector should contain the index for the corresponding PN orders, as indexed in posterior_dist_deltaphi[:, ..., :]), e.g. [[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]] for plotting all PN orders in the same plot (equivalent to nothing), or [[1], [2, 3, 4, 5], [6, 7, 8, 9, 10]] to obtain the LVK GWTC-3-like plot (figure 7).

In the config files,
"offset_x_axis_hierarchical_networks_upper_bounds": allows to offset the x-axis for the hierarchical upper bounds, so that they do not overlap if there are multiple networks [float, default 0.1]
"impose_y_axis_limits": allows to impose the y-axis limits for the plot, so that they are not automatically set based on the data [true/false]
"y_axis_limits": allows to set the y-axis limits for the plot, if impose_y_axis_limits is true [Nx2 vector, outer vector represent the N subplots, the inner vector represents the lower and upper limits, e.g. [[-1e-5, 1e-5],[-0.05, 0.05], [-1.0, 1.0]] ]
"violin_plots_bandwidth_std_multiplier": allows to set the bandwidth for the kernel density estimation used in the violin plots, as a multiplier of the standard deviation of the data [float, default 0.4]
"""
function plotDeltaPhiPosterior(plotTitle::AbstractString, posterior_dist_deltaphi::AbstractArray{Float64,3}; plot_conditioned_distribution::Bool = false, posterior_dist_deltaphi_conditioned::Union{Nothing, AbstractArray{Float64,3}} = nothing, subplots_pn_order_grouping::Union{Nothing, Vector{Vector{Int64}}} = nothing)
    
    set_common_plot_style()
    
    # Settings:
    plotHeight, plotWidth, plotDpi, padding = default_plot_dimensions() # Get the default plot dimensions

    # Define a list of different markers and colors to enhance readability
    markers, markersize, pointColors = get_markers_and_palette()

    # Proper LaTeX labels
    xlabel_str = "PN order" 
    ylabel_str = L"\delta \varphi_{\!i}"

    # Perform some checks and set defaults
    if plot_conditioned_distribution && posterior_dist_deltaphi_conditioned === nothing
        throw(ArgumentError("You set plot_conditioned_distribution to true, but posterior_dist_deltaphi_conditioned is nothing (should be a 3D array)."))
    end
    number_PN_orders_to_plot = size(posterior_dist_deltaphi, 2)
    if subplots_pn_order_grouping === nothing
        # If subplots_pn_order_grouping is nothing, then I will plot all PN orders in the same plot
        subplots_pn_order_grouping = [collect(1:number_PN_orders_to_plot)]
    end
    if any(x -> x < 1 || x > size(posterior_dist_deltaphi, 2), collect(Iterators.flatten(subplots_pn_order_grouping)))
        # Checks that all indices in subplots_pn_order_grouping are valid
        throw(ArgumentError("Invalid indices in subplots_pn_order_grouping."))
    end



    println("\n"*"#"^81)
    println("Plotting the posterior distribution for delta_phi.")

    # name_PN = collect(keys(pn_order_dic));
    PN_orders = configs["pn_waveforms"];
    name_PN = configs["pn_waveforms"];
    network_names = configs["network_list"];
    network_labels = labels_from_networks(network_names);
    PN_labels = labels_from_PN_orders(PN_orders);

    #Some additional checks
    if length(network_names) != size(posterior_dist_deltaphi, 1)
        throw(ArgumentError("The number of networks in the config file does not match the first dimension of posterior_dist_deltaphi."))
    end

    if length(PN_orders) != size(posterior_dist_deltaphi, 2)
        throw(ArgumentError("The number of PN orders in the config file does not match the second dimension of posterior_dist_deltaphi."))
    end
    if size(posterior_dist_deltaphi, 3) < 1
        throw(ArgumentError("The third dimension of posterior_dist_deltaphi should contain at least one sample."))
    end

    # Extend the arrays, eventually wrapping around
    markers = markers[mod1.(1:length(network_names), length(markers))]
    pointColors = pointColors[mod1.(1:length(network_names), length(pointColors))]

    # I will create a vector of subplots, as specified in subplots_pn_order_grouping, which then will be combined into a single plot
    subplots = Vector{Plots.Plot{Plots.GRBackend}}(undef, length(subplots_pn_order_grouping))

    small_margin_subplots = smaller_margins_subplots()  # Get the smaller margins for subplots

    # Loop over the subplots
    for (index_subplot, pn_order_indices) in enumerate(subplots_pn_order_grouping)
        # Create a new subplot for the current PN order grouping
        subplots[index_subplot] = plot(
            #xlabel=xlabel_str, 
            ylabel=(index_subplot > 1 ? "" : ylabel_str),
            #title=plotTitle,
            legend=(index_subplot == length(subplots_pn_order_grouping) ? :bottomright : false),
            xticks=(1:length(pn_order_indices), PN_labels[pn_order_indices]),
            #yscale=:log10, 
            size=(plotWidth * length(pn_order_indices) / number_PN_orders_to_plot, plotHeight), 
            dpi=plotDpi,
            xlims=(0.5 - configs["offset_x_axis_hierarchical_networks_upper_bounds"], length(pn_order_indices) + 0.5 + configs["offset_x_axis_hierarchical_networks_upper_bounds"]),
            grid=true,
            framestyle=:box,
            ylims=(configs["impose_y_axis_limits"] ? configs["y_axis_limits"][index_subplot] : :auto),  # Set y-axis limits if requested
            gridalpha=0.5, 
            gridcolor=:gray,  # Set grid lines to be transparent or gray
            yminorgrid=true, 
            minorgridalpha=0.3, 
            minorgridcolor=:gray,  # Enable minor grid lines
            xminorgrid=false,  # Disable minor grid lines for the x-axis,        
            left_margin = (index_subplot ==  1 ? small_margin_subplots.left_margin_with_ylabel : small_margin_subplots.left_margin),
            right_margin = small_margin_subplots.right_margin,
            top_margin = small_margin_subplots.top_margin,
            bottom_margin = small_margin_subplots.bottom_margin
        )

        # Add horizontal line at y = 0
        hline!(subplots[index_subplot], [0.], linecolor=:gray, linestyle=:dash, linewidth=2, label="")

        # Add top labels for each x-tick
        # for (pnindex, pno) in enumerate(pn_order_indices)
        #     annotate!(subplots[index_subplot], pnindex, maximum(posterior_dist_deltaphi[:, pno, :]) + 0.1, text(name_PN[pno], :center, 12))
        # end

        # Plot the posterior distributions for the current PN order grouping as violin plots, with some transparency and different color depending on the network
        for jj in 1:length(network_names)
            for (pnindex, pno) in enumerate(pn_order_indices)
                # plotting the violin plots for the delta_phi 'marginalized' distributions
                # It would be nice to have the width of the violin plots be proportional to the proability density function
                # but it seems that in Plot.jl there is no way to set the width of the violin plots
                violin!(
                    subplots[index_subplot], 
                    fill(pnindex + configs["offset_x_axis_hierarchical_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names)), length(posterior_dist_deltaphi[jj,pno,:])),
                    posterior_dist_deltaphi[jj,pno,:],
                    label=(pnindex == 1 ? network_labels[jj] : nothing),
                    #marker=markers[jj],  # Cycle through marker list
                    color=pointColors[jj], # Fill color
                    #markersize=6, 
                    #markerstrokewidth=1, 
                    #markerstrokecolor=:black,
                    orientation=:vertical,
                    alpha=0.45,          # Slight fill transparency for overlapping violin plots
                    linecolor=pointColors[jj],   # Outline color
                    linewidth=0,        # No outline thickness
                    width=0.65,         # Violin width... does not seem to work!
                    bandwidth=std(posterior_dist_deltaphi[jj,pno,:]) * configs["violin_plots_bandwidth_std_multiplier"] # Bandwidth for the kernel density estimation, to smooth out the (noisy) data
                )

                if plot_conditioned_distribution && posterior_dist_deltaphi_conditioned !== nothing
                    # Plot the conditioned distributions if requested
                    violin!(
                        subplots[index_subplot], 
                        fill(pnindex + configs["offset_x_axis_hierarchical_networks_upper_bounds"] * get_relative_x_offset_network(jj,length(network_names)), length(posterior_dist_deltaphi_conditioned[jj,pno,:])),
                        posterior_dist_deltaphi_conditioned[jj,pno,:],
                        label= nothing, # I set no labels at all here, just like in GWTC-3, you may explain in the caption that the conditioned distributions are plotted just as an outline
                        #marker=markers[jj],  # Cycle through marker list
                        color=:transparent, # Fill color
                        #markersize=10, 
                        #markerstrokewidth=10, 
                        #markerstrokecolor=:black,
                        orientation=:vertical,
                        linecolor=pointColors[jj],   # Outline color
                        linewidth= 10,        # Outline thickness -- does not seem to work!
                        width=0.65,         # Violin width, does not seem to work!
                        bandwidth=std(posterior_dist_deltaphi_conditioned[jj,pno,:]) * configs["violin_plots_bandwidth_std_multiplier"] # Bandwidth for the kernel density estimation, to smooth out the (noisy) data
                    )
                end
            end
       end

    end

    # Combine the subplots into a single plot
    final_plot = plot(subplots..., layout = @layout([grid(1, length(subplots_pn_order_grouping), widths = 0.999999999 .* [length(pn_indices) for pn_indices in subplots_pn_order_grouping] ./ number_PN_orders_to_plot )])) 
    #I will set the title in this final subplot, so that it is centered
    plot!(final_plot, title=plotTitle, xlabel=xlabel_str, size=(plotWidth, plotHeight), padding = padding, dpi=plotDpi) #, ylabel=ylabel_str, top_margin=2mm, bottom_margin=2mm, left_margin=2mm, right_margin=2mm)
    # Return the final plot
    return final_plot
end

