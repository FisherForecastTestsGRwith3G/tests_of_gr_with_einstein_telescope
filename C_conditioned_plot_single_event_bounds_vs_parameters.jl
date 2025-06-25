using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# Fast modification of the original script C, to produce a plot of the the trend of SNR and other parameters (as colored scatter plots) as a function the hierarchical upper limits, conditioned on the delta phi distribution.  

println("Running C_conditioned_plot_single_event_bounds_vs_parameters.jl script: this script is a modified version of the C_conditioned_plot.jl script, not completely refined. You should look a the source code and set the correct parameters in the script (not just in the config files), to obtain the wanted result.")

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

# Will be used to load the relevant parameters of the events
# gr_parameter_names = ["mc", "eta","chi1","chi2", "dL", "theta", "phi", "iota", "psi", "phiCoal", "tcoal", "lambda1", "lambda2"]


# For simplicity, also to not clutter too much the plot, I will iterate over a single detector network
list_indices_detector_network_to_use = [2]
# and for the same reason, I will iteratre over only the few specified PN orders (creating a separate plot for each of them)
list_indices_pn_order_to_use = [1]

# Create an empty array to store the results for the upper limits (mean and std for each of them... Eventually, if you have a single realization, the second parameter (std_dev) will be a NaN).
# upperLimits = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), 4) # 4 = mean, std_dev, mean_in_log_space, std_dev_in_log_space
# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]
upperLimitSingleEvents = zeros(length(list_of_networks), length(configs["pn_waveforms"]), n_events_used)

# Set up plots
set_common_plot_style()

# Settings:
plotHeight, plotWidth, plotDpi, padding = default_plot_dimensions() # Get the default plot dimensions

# Define a list of different markers and colors to enhance readability
markers, markersize, pointColors = get_markers_and_palette()

PN_orders = configs["pn_waveforms"];
network_names = configs["network_names"];
network_labels = labels_from_networks(network_names);
PN_labels = labels_from_PN_orders(PN_orders);

# process all the plots
 for (index_nn, nn) in enumerate(configs["network_list"])

    println("\n"*"#"^81)
    println("Processing Network = $(nn)\n")

    if !(index_nn in list_indices_detector_network_to_use)
        println("Skipping network $(nn), as it is not selected for the trend plot.")
        continue
    end

    # Loading the relevant parameters
    param_values = Dict()
    file_name = data_folder_name * "/catalog_w_deviations.h5"
    h5open(file_name, "r") do catalog_file
        param_group = catalog_file["parameter"]
        
        for param_name in keys(param_group)
            
            param_values[param_name] = read(param_group, param_name)
            #param_values_selected = param_values[global_index_total[nn]]
        end
    end


    #Produce the plot evenly spaced in log scale        
    numberOfEventsSingleRealization = n_events_realizations[index_realization]

    println("Processing realization = $(index_realization) out of $(number_of_points_trend_plot) with $(numberOfEventsSingleRealization) events")

    for (index_pno, pno) in enumerate(configs["pn_waveforms"])

        println("\n"*"#"^81)
        println("Processing PN = $(pno)\n")

        if !(index_pno in list_indices_pn_order_to_use)
            println("Skipping PN order $(pno), as it is not selected for the trend plot.")
            continue
        end
        
        # load the data
        folder = data_folder_name * nn * "/pn_" * pn_order_dic[pno][2]
        filename = folder * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # Load the SNR data
        snr_data = Dict()
        folder_name = data_folder_name * nn * "/" * pno_name *"/"
        file_name = folder_name * "snrs.h5"
        h5open(file_name, "r") do snr_file
            snr_data["values"] = read(snr_file, "values")  
            #snr_data["index"] = read(snr_file, "index")
        end

        vectorUpperLimits, upperLimitSingleEventsTemp, number_events_single_realization, numberOfEventsSingleRealizationUsed = obtain_conditioned_upper_bounds(
            configs["n_events"], 
            global_index_total[nn], 
            data["dphi0_k"], 
            data["delta_k"], 
            averageOverSeveralRealizations = false, 
            numberOfEventsSingleRealization = n_events_used, 
            n_events_and_numberOfEventsSingleRealization_refer_directly_to_observed_events = true,
            printEventsAsHorizontalLinesOrDensityPlot = true,
            use_all_n_events_for_single_event_sample_distribution = true,
            print_info_catalog_realization = false
        )

        upperLimitSingleEvents[index_nn, index_pno, 1:length(upperLimitSingleEventsTemp)] = upperLimitSingleEventsTemp



        # Produce the plot: a scatter plot with the single event upper bounds on the x axis, the SNR distribution on the y axis, and the param_values[selected_parameter] on the color axis.

        # Some quantities that may be useful to plot as the z axis
        
        total_mass_binary = @. param_values["mc"] / (param_values["η"])^(3. / 5.)  # Total mass of the binary system
        total_mass_binary = total_mass_binary[global_index_total[nn]]  # Select the global indices for the current network

        frequency_end_inspiral = @. 0.018 / (  param_values["mc"] / param_values["η"]^(3. /5.) ) / GMsun_over_c3
        frequency_end_inspiral = frequency_end_inspiral[global_index_total[nn]]

        # Quantities to be plotted
        x_axis_values = upperLimitSingleEventsTemp
        y_axis_values = snr_data["values"][global_index_total[nn]]
        z_axis_values = param_values["mc"][global_index_total[nn]]
        # z_axis_values = total_mass_binary
        
        # Proper LaTeX labels
        xlabel_str =L"|\delta" * LaTeXString(PN_labels[pno]) * L"|", 
        ylabel_str = "SNR" 
        zlabel_str = L"\mathcal{M}_c"  # or L"M_{tot}" for total mass

        # Check that the length of the x, y, and z axis values are the same
        if length(x_axis_values) != length(y_axis_values) || length(x_axis_values) != length(z_axis_values)
            error("The lengths of x, y, and z axis values must be the same. Got lengths: x=$(length(x_axis_values)), y=$(length(y_axis_values)), z=$(length(z_axis_values))")
        end

        # Now I produce the scatter plot, in log for all three axis
        final_plot = scatter(
            x_axis_values, 
            y_axis_values, 
            z_axis_values,
            xlabel=xlabel_str,
            ylabel=ylabel_str, 
            zlabel=zlabel_str,  # or L"M_{tot}"
            title="",
            #color=:blues,  # Use a color gradient
            markersize=5,
            marker=:circle,
            xscale=:log10, 
            yscale=:log10, 
            zscale=:log10,
            #legend=:topright
            size=(plotWidth, plotHeight), 
            # padding = padding, 
            # dpi=plotDpi,
            # framestyle = :box,    
            # minorgridcolor=:gray,  # Enable minor grid lines
            # xminorgrid=true,  # Disable minor grid lines for the x-axis,     
            # grid = true,
            # gridalpha=0.4, 
            # gridcolor=:gray,  # Set grid lines to be transparent or gray
            # yminorgrid=true, 
            # minorgridalpha=0.15
        )

        # Save the plot
        output_folder_name_plot = output_folder_name * nn * "/parameters_vs_conditioned_delta_phi_upper_limits/"
        mkpath(output_folder_name_plot)
        plot_filename = output_folder_name_plot * "plot_parameters_vs_conditioned_delta_phi_upper_limits_" * simulation_tag * "_" * nn * "_" * pno_name * ".pdf"
        println("Saving the conditioned upper limits plot to file: ", plot_filename)
        savefig(final_plot, plot_filename)

    end
end