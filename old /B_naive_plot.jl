using HDF5
using Plots
import Contour
using Trapz
using Plots.PlotMeasures

include(".paths.jl")
default(fontfamily="Computer Modern", labelfontsize=21, titlefontsize=19, guidefontsize=10, tickfontsize=16, legendfontsize=12)
using LaTeXStrings
using Statistics

################################################################################
## Specify simulation specs in this part of the script
# ||||||||
# vvvvvvvv

PhD = "Andrea"
path_catalog, path_output = whoIsThere(PhD)

# specify where the data is stored
simulation_tag = "BGR_HM"

# network specs
network_names = #["ETS"]
                ["network_0_15km", "network_45_15km", "ETS"]
                
# specify the pn-orders we want to analze
pn_orders = #["0", "0.5", "1", "1.5"] 
            ["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"]   

# snr threshold
snr_thresh = 12.

# secify where to save the plots
figure_dir = path_output*"output/"*simulation_tag*"/plots/"

### Restrict_catalog 
n_events = 10000  # number of events from the catalog used     

# ^^^^^^^^
# ||||||||
## Specify simulation specs in this part of the script
################################################################################

function plot_(mean_cumulative_error)

    # Use LaTeX-like fonts

    # Proper LaTeX labels
    xlabel_str = L"\mathrm{PN~Order}"  # Use \mathrm for proper LaTeX rendering
    ylabel_str = L"\mathrm{Cumulative~Error}"

    # Define a list of different markers to enhance readability
    markers = [:circle, :square, :diamond, :utriangle, :dtriangle, :hexagon]

    # Calculate mean and standard deviation for each ii and jj
    # mean_cumulative_error = zeros(length(pn_orders), length(network_names))
    # std_cumulative_error = zeros(length(pn_orders), length(network_names))

    plot_padding_top= 3;
    plot_padding_bottom = 4;

    # Find the minimum and maximum values for y-axis limits, excluding the first point
    min_y = minimum(mean_cumulative_error[2:end, :]) / plot_padding_bottom  # Divide by 10 for some padding
    max_y = maximum(mean_cumulative_error[2:end, :]) * plot_padding_top  # Multiply by 10 for some padding

    # Find the minimum and maximum values for y-axis limits, excluding the first point
    min_y_sub = minimum(mean_cumulative_error[1, :]) / plot_padding_bottom  # Divide by 10 for some padding
    max_y_sub = maximum(mean_cumulative_error[1, :]) * plot_padding_top  # Multiply by 10 for some padding

    # println("min_y ", min_y)
    # println("max_y ", max_y)
    # println("min_y_sub ", min_y_sub)
    # println("max_y_sub ", max_y_sub)
    # println(size(mean_cumulative_error))

    # Initialize the main plot with white background and high resolution
    cc_main = plot(
        xlabel=xlabel_str, title="Cumulative Error on PN Order - PhenomHM ET ",
        legend=:bottomright, xticks=(2:10, name_PN[2:end]), 
        yscale=:log10, size=(900, 600), dpi=300,
        xlims=(1.5, 10.5),
        grid=true, framestyle=:box,
        yticks=[10.0^i for i in floor(Int, log10(min_y)):ceil(Int, log10(max_y))],  # Set y-axis ticks for each 10^N value within the range
        ylims=(min_y, max_y),  # Automatically set the y-axis range
        gridalpha=0.5, gridcolor=:gray,  # Set grid lines to be transparent or gray
        yminorgrid=true, minorgridalpha=0.3, minorgridcolor=:gray,  # Enable minor grid lines
        xminorgrid=false  # Disable minor grid lines for the x-axis
    )
    display(cc_main)

    # Initialize the subplot for the first point
    cc_sub = plot(
        xlabel=xlabel_str, ylabel=ylabel_str,
        legend=false, xticks=([1], [name_PN[1]]), 
        xlims=(0.5, 1.5),
        yscale=:log10, size=(300, 600), dpi=300,
        grid=true, framestyle=:box,
        yticks=[10.0^i for i in floor(Int, log10(min_y_sub)) : ceil(Int, log10(max_y_sub))],  # Set y-axis ticks for the first point
        ylims=(min_y_sub, max_y_sub),  # Automatically set the y-axis range
        gridalpha=0.5, gridcolor=:gray,  # Set grid lines to be transparent or gray
        yminorgrid=true, minorgridalpha=0.3, minorgridcolor=:gray,  # Enable minor grid lines
        xminorgrid=false  # Disable minor grid lines for the x-axis
    )

    display(cc_main)

    # # Plot the first point in the subplot
    for jj in 1:length(network_names)
        scatter!(
            cc_sub, [1], [mean_cumulative_error[1, jj]],
            label=network_names[jj],
            marker=markers[mod1(jj, length(markers))],  # Cycle through marker list
            markersize=6, markerstrokewidth=1, markerstrokecolor=:black,
            alpha=0.9  # Slight transparency for overlapping points
        )
    end

    # Plot the remaining points in the main plot
    for jj in 1:length(network_names)
        scatter!(
            cc_main, 2:length(pn_orders), mean_cumulative_error[2:end, jj],
            label=network_names[jj],
            marker=markers[mod1(jj, length(markers))],  # Cycle through marker list
            markersize=6, markerstrokewidth=1, markerstrokecolor=:black,
            alpha=0.9  # Slight transparency for overlapping points
        )
        
    #     # Plot individual errors as small transparent orange horizontal stripes
    #     for ii in 1:1
    #         vecc = PN_errors[ii, jj, 1:NUMBER_OF_EVENTS_TO_BE_USED_SINGLE_CATALOG]
    #         scatter!(
    #             cc_sub, fill(ii, length(vecc)), vecc,
    #             marker=:hline, markersize=15, markerstrokewidth=2, alpha=0.5, color=:orange, label=""
    #         )
    #     end

    #     for ii in 2:length(pn_orders)
    #         vecc = PN_errors[ii, jj, 1:NUMBER_OF_EVENTS_TO_BE_USED_SINGLE_CATALOG]
    #         scatter!(
    #             cc_main, fill(ii, length(vecc)), vecc,
    #             marker=:hline, markersize=15, markerstrokewidth=2, alpha=0.5, color=:orange, label=""
    #         )
    #     end
    end

    # Add a dummy plot for the legend entry in the main plot
    # scatter!(
    #     cc_main, [NaN], [NaN],
    #     marker=:hline, markersize=15, markerstrokewidth=2, alpha=0.5, color=:orange, label="Errors from single events"
    # )

    # # Plot the GWTC-3 results in the main plot
    # scatter!(
    #     cc_main, 2:length(pn_orders), LVK_GWTC3_results[2:end],
    #     label="LVK GWTC-3", marker=:diamond, markersize=6, markerstrokewidth=1, markerstrokecolor=:black, color=:darkblue
    # )

    # # Plot the first GWTC-3 result in the subplot
    # scatter!(
    #     cc_sub, [1], [LVK_GWTC3_results[1]],
    #     label="LVK GWTC-3", marker=:diamond, markersize=6, markerstrokewidth=1, markerstrokecolor=:black, color=:darkblue
    # )

    # Combine the main plot and the subplot
    final_plot = plot(cc_sub, cc_main, layout = @layout [a{0.15w} b{0.85w}])#, top_margin=2mm, bottom_margin=2mm, left_margin=2mm, right_margin=2mm)

    # Save the combined plot to a file on the remote server

    # Display the final plot
    return final_plot
end




#######################################################################################################################################
output_folder_name =  path_output*"output/"*simulation_tag*"/data/"
name_PN = ["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"];

pn_order_dic = Dict(
    "-1"       => (-1.0    ,"minus_one"),
    "0"        => (0.0     ,"zero"), 
    "0.5"      => (0.5     ,"half"),
    "1"        => (1.0     ,"one"), 
    "1.5"      => (1.5     ,"one_half"),
    "2"        => (2.0     ,"two"),
    "log(2.5)" => (log(2.5),"log_two_half"),
    "3"        => (3.0     ,"three"), 
    "log(3.)"  => (log(3.) ,"log_three"),
    "3.5"      => (3.5     ,"three_half")
);

### get global indices 
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()
for nn in network_names

    file_name = output_folder_name * nn *"/global_indices.h5"
    h5open(file_name, "r") do gi_file
        global_index_fisher[nn] = read(gi_file, "fisher")  
        global_index_snr[nn] = read(gi_file, "snr")
        global_index_total[nn] = read(gi_file, "total")
    end
end


# cumulative error on PNorder
cumulative_error = zeros(length(pn_orders), length(network_names))


mkpath(figure_dir*"/naive_plot/")

### analyze catalog for each selected pn order
event_ks = collect(1:n_events)
global ii = 0
for nn in network_names
    global ii += 1
    idx_g_total = global_index_total[nn][1:n_events]
    idx_g_fisher = global_index_fisher[nn][1:n_events]

    global jj = 0
    for pno in pn_orders
        global jj += 1
        pno_name = "pn_" *pn_order_dic[pno][2]
        folder_name = output_folder_name * nn * "/" * pno_name *"/"

        #load snr
        snr_data = Dict()
        file_name = folder_name * "snrs.h5"
        h5open(file_name, "r") do snr_file
            snr_data["values"] = read(snr_file, "values")  
            snr_data["index"] = read(snr_file, "index")
        end

        #load estimated parameter and deviations 
        data = Dict()
        file_name = folder_name *"single_event_measurements.h5"
        h5open(file_name, "r") do data_file
            data["dphit_k"] = read(data_file, "dphit_k")  
            data["dphi0_k"] = read(data_file, "dphi0_k")
            data["delta_k"] = read(data_file, "delta_k")
        end

        # plot snr 
        snr_values = snr_data["values"][1:n_events]
        
        # plot fisher error summary
        fisher_errors = data["delta_k"][1:n_events]


        fisher_errors_non_Nan = fisher_errors[.!isnan.(fisher_errors)]

        cumulative_error[jj,ii] = sum(fisher_errors_non_Nan.^(-2))^-0.5
    end
end
plot_HM = plot_(cumulative_error)

println("Plotting cumulative error on PN order in ", string(figure_dir*"naive_plot/"*"cumulative_error.png"))
savefig(plot_HM, figure_dir*"naive_plot/"*"cumulative_error.pdf")


