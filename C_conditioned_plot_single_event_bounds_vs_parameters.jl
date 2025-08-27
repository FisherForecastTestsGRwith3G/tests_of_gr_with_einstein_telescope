using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing
using Printf

# Needed for the get_dL function
using Pkg
parentdir = dirname(pwd())
Pkg.activate(string(parentdir)*"/GW.jl")
using GW: get_dL

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
list_indices_detector_network_to_use = [1,2,3]
# and for the same reason, I will iteratre over only the few specified PN orders (creating a separate plot for each of them)
# list_indices_pn_order_to_use = [1]
list_indices_pn_order_to_use = collect(1:length(configs["pn_waveforms"])) # Analyzes all PN orders

# Create an empty array to store the results for the upper limits (mean and std for each of them... Eventually, if you have a single realization, the second parameter (std_dev) will be a NaN).
# upperLimits = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), 4) # 4 = mean, std_dev, mean_in_log_space, std_dev_in_log_space
# Create an empty array to store the single events upper limits, if required.
n_events_used = configs["n_events"]
upperLimitSingleEvents = zeros(length(configs["network_list"]), length(configs["pn_waveforms"]), n_events_used)

# Set up plots
set_common_plot_style()

# Settings:
plotHeight, plotWidth, plotDpi, padding = default_plot_dimensions() # Get the default plot dimensions

# Define a list of different markers and colors to enhance readability
markers, markersize, pointColors = get_markers_and_palette()

PN_orders = configs["pn_waveforms"];
network_names = configs["network_list"];
network_labels = labels_from_networks(network_names);
PN_labels = labels_from_PN_orders(PN_orders);

# not used in the end, since Plots.jl does not support log colorbars nor costum colorbar ticks
function auto_log_colorbar_ticks(z_axis_values; min_ticks=4, max_ticks=8)
    # Get min and max in linear space
    zmin = minimum(z_axis_values)
    zmax = maximum(z_axis_values)
    # Convert back to linear
    zmin_lin = 10.0^zmin
    zmax_lin = 10.0^zmax

    # Find the nearest powers of 10
    pow_min = floor(Int, log10(zmin_lin))
    pow_max = ceil(Int, log10(zmax_lin))

    # Try different intervals
    for interval in log10.((1000, 500, 200, 100, 50, 20, 10,5,2, 1, 0.5, 0.2, 0.1))
        ticks = collect(pow_min:interval:pow_max)
        if min_ticks <= length(ticks) <= max_ticks
            tick_positions = ticks
            tick_labels = [@sprintf("%.0f", 10.0^t) for t in ticks]
            return tick_positions, tick_labels
        end
    end

    # Fallback: just use 5 ticks
    ticks = range(pow_min, pow_max; length=5)
    tick_positions = collect(ticks)
    tick_labels = [@sprintf("%.0f", 10.0^t) for t in ticks]
    return tick_positions, tick_labels
end

function get_z(dL)
    # use bisection method to find z
    z=1
    zmin = 0.
    zmax = 30.
    dL = dL*1e3 # convert to Mpc
    iteration = 0
    while abs.(dL - get_dL(z)[1]) .> 1e-3
        iteration += 1
        if iteration > 100
            println("Could not find z")
            break
        end
        z = 0.5(zmin + zmax)
        if dL .> get_dL(z)[1]
            zmin = z
        else
            zmax = z
        end
    end
    return z
end

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
    file_name = user_configs["path_output"]*simulation_tag* "/catalog_w_deviations.h5"
    h5open(file_name, "r") do catalog_file
        param_group = catalog_file["parameter"]
        
        for param_name in keys(param_group)
            
            param_values[param_name] = read(param_group, param_name)
            #param_values_selected = param_values[global_index_total[nn]]
        end
    end


    #Produce the plot evenly spaced in log scale   

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

        # Evaluate some derived quantities
        param_values["z"] = get_z.(param_values["dL"]) 
        param_values["total_mass_binary"] = @. param_values["mc"] / (param_values["eta"])^(3. / 5.)  # Total mass of the binary system (in detector frame)
        GMsun_over_c3 = 4.925491025543575903411922162094833998e-6 # seconds
        param_values["tfrequency_end_inspiral"] = @. 0.018 / (  param_values["mc"] / param_values["eta"]^(3. /5.) ) / GMsun_over_c3 # Frequency at the end of the inspiral (in Hz), in detector frame
        param_values["mc_source_frame"] = @. param_values["mc"] / (1 + param_values["z"])  # Chirp mass in source frame
        param_values["total_mass_binary_source_frame"] = @. param_values["total_mass_binary"] / (1 + param_values["z"])  # Chirp mass in source frame
        param_values["q"] = @. -((-1. + sqrt(1. - 4. * (param_values["eta"])) + 2. * (param_values["eta"]) ) / (2. * (param_values["eta"]) )) # mass ratio q = m2 / m1 
        param_values["m1"] = @. param_values["mc"] * ((1. + param_values["q"] ) / param_values["q"]^3)^(1. / 5.) # detector frame mass of the primary component
        param_values["m2"] = @. param_values["m1"] * param_values["q"] # detector frame mass of the secondary component
        param_values["m1_source_frame"] = @. param_values["m1"] / (1 + param_values["z"]) # source frame mass of the primary component
        param_values["m2_source_frame"] = @. param_values["m2"] / (1 + param_values["z"]) # source frame mass of the secondary component
        
        # Load the SNR data
        snr_data = Dict()
        pno_name = "pn_"*pn_order_dic[pno][2]
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
        
        # Whether to overlay the scatter plot with markers for ill-conditioned events
        # This does not make sense in the current setup below (so plot_ill_conditioned_events = false right now), but would for example make sense when plotting iota and mass ratio on x and y axis, with the PN bound coefficients on the z axis
        plot_ill_conditioned_events = false 

        # Quantities to be plotted
        x_axis_values = upperLimitSingleEventsTemp
        y_axis_values = snr_data["values"]
        z_axis_values = param_values["mc"]
        # z_axis_values = total_mass_binary

        # Another set of parameters that could have been plotted (and would require disabling the SNR dashed line below)
        # x_axis_values = param_values["iota"]
        # y_axis_values = 1. ./ param_values["q"]
        # z_axis_values = upperLimitSingleEventsTemp

        if plot_ill_conditioned_events
            # Saves also the values for the events that should have been used, but have not been used, due to ill-conditioned Fisher matrix
            # So the events with global_index_snr[nn] set to true but global_index_fisher[nn] set to false
            indices_ill_conditioned = global_index_snr[nn] .& .!global_index_fisher[nn]
            x_ill_conditioned = x_axis_values[indices_ill_conditioned]
            y_ill_conditioned = y_axis_values[indices_ill_conditioned]
            #z_ill_conditioned = z_axis_values[indices_ill_conditioned]
        end

        # Obtain correct array elements
        # x_axis_values = x_axis_values[global_index_total[nn]]
        y_axis_values = y_axis_values[global_index_total[nn]]
        z_axis_values = z_axis_values[global_index_total[nn]]
        
        # Proper LaTeX labels
        xlabel_str =L"|\delta" * LaTeXString(PN_labels[index_pno]) * L"|"
        ylabel_str = "SNR" 
        zlabel_str = L"\mathcal{M}_c"  # or L"M_{tot}" for total mass

        # Other set of possible parameters for example
        # xlabel_str =L"\iota"
        # ylabel_str = L"1/q" 
        # zlabel_str = L"|\delta" * LaTeXString(PN_labels[index_pno]) * L"|"  # or L"M_{tot}" for total mass

        # Check that the length of the x, y, and z axis values are the same
        # The code above and below must be accordingly modified, since some arrays just contain the observed events (like upperLimitSingleEventsTemp), some other contain all the events (like snr_data["values"] and param_values["mc"])
        if length(x_axis_values) != length(y_axis_values) || length(x_axis_values) != length(z_axis_values)
            error("The lengths of x, y, and z axis values must be the same. Got lengths: x=$(length(x_axis_values)), y=$(length(y_axis_values)), z=$(length(z_axis_values))")
        end

        use_log_scale_z_axis = true

        if use_log_scale_z_axis
            println("Using log scale for the z axis. Plots.jl does not correctly support log scale, nor arbitrary colorbar ticks, so there is no concrete easy way to display log data on the z axis. I will just plot directly log_10() of the quantity")
            z_axis_values = log10.(z_axis_values)
            zlabel_str = L"\log_{10}("*LaTeXString(zlabel_str)*L")"
            tick_positions, tick_labels = auto_log_colorbar_ticks(z_axis_values)
        end

        # Print log 10 of minimum and maximum values of the z axis, if use_log_scale_z_axis is true, otherwise print the minimum and maximum values
        if use_log_scale_z_axis
            println("Minimum value of z axis (log10): ", minimum(z_axis_values))
            println("Maximum value of z axis (log10): ", maximum(z_axis_values))
        else
            println("Minimum value of z axis: ", minimum(z_axis_values))
            println("Maximum value of z axis: ", maximum(z_axis_values))
        end

        # Now I produce the scatter plot, in log for all three axis
        final_plot = scatter(
            x_axis_values, 
            y_axis_values, 
            marker_z = z_axis_values,
            # xlims = (xmin, xmax),
            # ylims = (ymin, ymax),
            # clims = (cmin, cmax), 
            color = :viridis,
            xlabel=xlabel_str,
            ylabel=ylabel_str, 
            colorbar = true,
            colorbar_title=zlabel_str,  # or L"M_{tot}"
            title="",
            label = "Single events",
            xscale=:log10, 
            yscale=:log10, 
            #zscale=:log10,
            #colorbar_scale=:log10, # bug in Plots.jl, so I will not use it...
            #zscale=:log10,
            legend= :topright,
            size=(plotWidth, plotHeight),
            right_margin = 5mm, # Increase right_margin to avoid cut off labels 
            padding = padding*1.2,  # Increase padding to avoid cut off labels
            dpi=plotDpi*0.6, # Reduce DPI to decrease file size
            framestyle = :box,    
            minorgridcolor=:gray,  # Enable minor grid lines
            xminorgrid=true,  # Disable minor grid lines for the x-axis,     
            grid = true,
            gridalpha=0.6, 
            gridcolor=:gray,  # Set grid lines to be transparent or gray
            yminorgrid=true, 
            minorgridalpha=0.2,
            markersize=2,
            marker=:circle,
            alpha=0.45,  # Set transparency for the markers
            markerstrokewidth=0,
            #colorbar_ticks = (tick_positions, tick_labels) # Not implemented in Plots.jl, so I will not use it...
        )

        if plot_ill_conditioned_events
            # Overlay the scatter plot with purple "X" markers for ill-conditioned events
            # The code above must be changed to have two parameters on the x and y which are evaluated before the Fisher matrix
            # For example it makes sense to set plot_ill_conditioned_events = true when iota and mass ratio are on x and y axis, with the PN bound coefficients on the z axis
            scatter!(
                x_ill_conditioned,
                y_ill_conditioned,
                color = :darkorange,
                marker = :x,
                markersize = 2,
                alpha=0.45,  # Set transparency for the markers
                label = "Events with ill-conditioned FIM"
            )
        end

        # overlay an horizontal dashed gray line at configs["snr_thresh"], which should be SNR = 12
        # Keep only if applicable!
        hline!(final_plot, [configs["snr_thresh"]], linestyle=:dash, color=:gray, alpha = 1.0, label="SNR threshold = $(Int(configs["snr_thresh"]))")

        # Save the plot
        output_folder_name_plot = output_folder_name * "parameters_vs_conditioned_delta_phi_upper_limits/" * nn * "/"
        mkpath(output_folder_name_plot)
        # I do not use pdfs but png instead, since there are about 150.000 points in each plot, and therefore the pdfs would be slow to load.
        plot_filename = output_folder_name_plot * "plot_parameters_vs_conditioned_delta_phi_upper_limits_" * simulation_tag * "_" * nn * "_" * pno_name * ".png"
        println("Saving the conditioned upper limits plot to file: ", plot_filename)
        savefig(final_plot, plot_filename)

    end
end