using HDF5
using Plots
using Trapz
using LaTeXStrings
using Random 
using Base.Threads
using Plots.Measures

# include required scripts 
include("_hierarchical_dist.jl")
include("_hierarchical_plotting_utils.jl")
include("_parse_config.jl")
include("_plot_style.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

# TODO: Write a better config reader 
configs, simulation_tag = readConfigForE(config_file_name)

# set up folder names
catalog_w_deviations_folder = user_configs["path_output"]*simulation_tag*"/"
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

@threads for nn in configs["network_list"]

    contours = Dict()
    
    for idx_noise in 0:configs["n_noise_realizations"]     
        
        noise_seed = configs["noise_seed"] + idx_noise -1  
        rng =  Xoshiro(noise_seed) 

        if idx_noise == 0
            noise_seed = "stored_data"
        end
         
        ### process all the contours
        
        contours[noise_seed] = Dict()

        # initiate the limits of the plots
        mu_limit_overlay = [0., 0.]
        sig_limit_overlay = [0., 0.]
        mu_limit_overlay_inset = [0., 0.]
        sig_limit_overlay_inset = [0., 0.]

        for pno in configs["pn_waveforms"]

            println("\n"*"#"^81)
            println("Processing PN = $(pno), Noise seed = $(noise_seed)\n")
            
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
            n_events_used = configs["n_events"]
            if n_events_used > sum(global_index_total[nn])
                println("Warning: There are at most $(sum(global_index_total[nn]))-datapoints available")
                n_events_used = sum(global_index_total[nn])
            else
                println("Events used: $(n_events_used)")
            end
            
            dphit_k = data["dphit_k"][global_index_total[nn]]
            dphit_k = dphit_k[1:n_events_used] 
            dphi0_k = data["dphi0_k"][global_index_total[nn]]
            dphi0_k = dphi0_k[1:n_events_used] 
            delta_k = data["delta_k"][global_index_total[nn]]
            delta_k = delta_k[1:n_events_used]
            if idx_noise > 0 # generate a new noise realization
                dphi0_k = dphit_k .+ delta_k .* randn(rng, n_events_used)
            end

            ### calculate the hyper-parameter distribution
            # first estimate where to place it
            center_mu = sum(dphi0_k) / n_events_used
            center_sig = max(0, sqrt.(sum((center_mu .- dphi0_k).^2) ./  n_events_used))
            spread = sqrt.(1.0 ./ sum(1 ./ (center_sig.^2 .+ delta_k.^2) ))
            println("center_sig: ", center_sig)
            println("spread: ", spread)
            println("center_mu: ", center_mu)
            #spread = sqrt.(1.0 ./ sum(1 ./ delta_k.^2) )


            k_spread = configs["k_spread"]
            mu_limit = (center_mu - k_spread*spread, center_mu + k_spread*spread)
            if !isnothing(configs["mu_lims"][pno]) # overwrite mu_lims
                mu_limit = (configs["mu_lims"][pno][1], configs["mu_lims"][pno][2])
                
            end 
            mu_limit_overlay[1] = min(mu_limit[1], mu_limit_overlay[1])
            mu_limit_overlay[2] = max(mu_limit[2], mu_limit_overlay[2])
            if pno in configs["inset_pns"]
                mu_limit_overlay_inset[1] = min(mu_limit[1], mu_limit_overlay_inset[1])
                mu_limit_overlay_inset[2] = max(mu_limit[2], mu_limit_overlay_inset[2])
            end
            println("Mu-limits: $(mu_limit)")
            
            sig_limit = (0, center_sig + k_spread*spread)
            if !isnothing(configs["sigma_lims"][pno]) # overwrite sig_lims
                sig_limit = (configs["sigma_lims"][pno][1], configs["sigma_lims"][pno][2])
                
            end 
            sig_limit_overlay[1] = min(sig_limit[1], sig_limit_overlay[1])
            sig_limit_overlay[2] = max(sig_limit[2], sig_limit_overlay[2])
            if pno in configs["inset_pns"]
                sig_limit_overlay_inset[1] = 0.
                sig_limit_overlay_inset[2] = max(sig_limit[2], sig_limit_overlay_inset[2])
            end
            println("Sigma-limits: $(sig_limit)")

            # calculate the ranges 
            mu_values = collect(LinRange(mu_limit[1], mu_limit[2], configs["n_points"]))
            sig_values = collect(LinRange(sig_limit[1], sig_limit[2], configs["n_points"]))
            p_mu_sig, p_sig, p_mu, n_tot, nn_marg = hyperparamDistTIGER(mu_values, sig_values, dphi0_k, delta_k)
            
            println("Consistency check for calculated distribution:")
            println("Difference in normalization = $(n_tot-nn_marg)")
            
            # calculate contour
            contour_lines_x, contour_lines_y = calculate90CIContour(mu_values, sig_values, p_mu_sig) 
            contours[noise_seed][pno] = [contour_lines_x, contour_lines_y]
        end

        # create the overlayed plot 
        set_common_plot_style()

        overlayed_contour_plot = plot([0,0], sig_limit_overlay, linestyle = :dash, linecolor=:gray, label=false, legend=:topleft, left_margin = 3mm, bottom_margin = 3mm, right_margin = 3mm, top_margin = 3mm, size = (800,600))
        plot!(overlayed_contour_plot, [0,0], sig_limit_overlay, linestyle = :dash, linecolor=:gray, label=false, inset=bbox(0.05, 0.05, 0.35, 0.35, :right), subplot=2, legend=false,)
        zoomed_in_contour_plot = overlayed_contour_plot[2]
        for pno in configs["pn_waveforms"]
            n_countours = length(contours[noise_seed][pno][1])
            label = labels_from_pn_orders[pno]
            color = color_from_pn_orders[pno]
            

            if pno == "-1"
                label = label * L"\times(500)"
            end
            for idx_c in 1:n_countours
                x_data = contours[noise_seed][pno][1][idx_c]
                y_data = contours[noise_seed][pno][2][idx_c]
                if pno == "-1"
                    x_data = x_data*500
                    y_data = y_data*500
                end
                plot!(
                    overlayed_contour_plot[1], 
                    x_data, 
                    y_data, 
                    label = label,
                    color = color,
                    linewidth=2)
            end
        end

        xlabel!(overlayed_contour_plot[1], "μ")
        ylabel!(overlayed_contour_plot[1],"σ")
        xlims!(overlayed_contour_plot[1], mu_limit_overlay...)
        ylims!(overlayed_contour_plot[1], sig_limit_overlay...)
        plot!(overlayed_contour_plot[1], right_margin = 3mm)
        plot!(overlayed_contour_plot[1], alpha = 1.)

        # create zoomed in overlayed plot
        #zoomed_in_contour_plot = plot([0,0], sig_limit_overlay_inset, linestyle = :dash, linecolor=:gray, label=false)
        for pno in configs["pn_waveforms"]
            n_countours = length(contours[noise_seed][pno][1])
            label = labels_from_pn_orders[pno]
            color = color_from_pn_orders[pno]
            

            if pno == "-1"
                label = label * L"\times(500)"
            end
            for idx_c in 1:n_countours
                x_data = contours[noise_seed][pno][1][idx_c]
                y_data = contours[noise_seed][pno][2][idx_c]
                if pno == "-1"
                    x_data = x_data*500
                    y_data = y_data*500
                end
                plot!(
                    zoomed_in_contour_plot, 
                    x_data, 
                    y_data, 
                    color = color,
                    linewidth=2,
                    xticks = [
                        mu_limit_overlay_inset[1], 
                        0.,
                        mu_limit_overlay_inset[2]])
            end
        end

        xlims!(zoomed_in_contour_plot, mu_limit_overlay_inset...)
        ylims!(zoomed_in_contour_plot, sig_limit_overlay_inset...)
        plot!(zoomed_in_contour_plot, alpha = 1.)

        # save plot
        mkpath(output_folder_name* nn * "/all_orders/")
        fig_file_name = output_folder_name * nn * "/all_orders/hyperparamdist_90CI_$(noise_seed).pdf"
        println("\nSaving 90CI plot in: $(fig_file_name)")
        savefig(overlayed_contour_plot, fig_file_name)
    end

end

