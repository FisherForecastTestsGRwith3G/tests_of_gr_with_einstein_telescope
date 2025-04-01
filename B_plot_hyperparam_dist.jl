using HDF5
using Plots
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing

# include required scripts 
include("_hierarchical_dist.jl")
include("_hierarchical_plotting_utils.jl")
include(".paths.jl")
include("_setup_MCMC.jl")

################################################################################
## Specify simulation specs in this part of the script
# ||||||||
# vvvvvvvv

include("_simulation_settings.jl")

# Below here you may ovveride the simulation settings specified in the _simulation_settings.jl file
# Yet for consistency across different scripts it is recommended to set the simulation settings in the _simulation_settings.jl file

# specify the pn-orders we want to analze
pn_orders =  ["0"]
           # ["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"]   

#MCMC = true

output_folder_name = output_folder_name*"data/"

### plot specifiers
n_events = 30   # number of events from the catalog used 

# quick and dirty limit option
n_points = 1000 # Determines gridpoints for visualizing the distribution
                # Higher value improves estimate of evidence and percentiles
                # but also increases computing time.
k_spread = 30   # multipies the estimated spread of the distribution
                # for automatic setting of the plotting limits. 
                # If limits dont make sense, increase this value as 
                # a first quick fix.

# set limits manually
mu_lims_dict = Dict(
    "-1"       => nothing,  #<- if set to nothing, limit will be
    "0"        => (-0.05, 0.05),  # tried to be infered outmatically
    "0.5"      => (-0.05, 0.05),  # Unsing k_spread variable
    "1"        => (-0.05, 0.05),  # Alternatively, the limits can 
    "1.5"      => (-0.05, 0.05),  # be specifeid as tuples i.e.
    "2"        => nothing,  # (mu_lim_low, mu_lim_up)
    "log(2.5)" => nothing,
    "3"        => nothing,
    "log(3.)"  => nothing,
    "3.5"      => nothing
    )   

sig_lims_dict = Dict(
    "-1"       => nothing,
    "0"        => (0.0, 0.1),
    "0.5"      => (0.0, 0.1),
    "1"        => (0.0, 0.1),
    "1.5"      => (0.0, 0.1),
    "2"        => nothing,
    "log(2.5)" => nothing,
    "3"        => nothing,
    "log(3.)"  => nothing,
    "3.5"      => nothing
    );

# injected values
val_injected_dict = Dict(
    "-1"       => (0.,0.0025),  #<- plots lines for injected values
    "0"        => (0.,0.0025),  # if set to noting, no lines are plotted
    "0.5"      => (0.,0.0025),  # Otherwise, specify as tuple i.e.
    "1"        => (0.,0.0025),  # (mu_injected, sig_injected)  
    "1.5"      => (0.,0.0025),
    "2"        => (0.,0.0025),
    "log(2.5)" => (0.,0.0025),
    "3"        => (0.,0.0025),
    "log(3.)"  => (0.,0.0025),
    "3.5"      => (0.,0.0025)
    )             

# ^^^^^^^^
# ||||||||
## Specify simulation specs in this part of the script
################################################################################


# get global index
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

# process all the plots
 for nn in network_names

    for pno in pn_orders

        println("\n"*"#"^81)
        println("Processing PN = $(pno)\n")
        
        # load the data
        folder = output_folder_name * nn * "/pn_" * pn_order_dic[pno][2]
        filename = folder * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # extract only valid data-points via global index
        n_events_used = n_events
        if n_events > sum(global_index_total[nn])
            println("Warning: There are at most $(sum(global_index_total[nn]))-datapoints available")
            n_events_used = sum(global_index_total[nn])
        end
        dphi0_k = data["dphi0_k"][global_index_total[nn]]
        dphi0_k = dphi0_k[1:n_events_used] 
        delta_k = data["delta_k"][global_index_total[nn]]
        delta_k = delta_k[1:n_events_used]

        ### calculate the hyper-parameter distribution
        # first estimate where to place it
        center_mu = sum(dphi0_k) / n_events_used
        center_sig = max(0, sqrt.(sum(center_mu .- dphi0_k).^2 ./  n_events_used))
        spread = sqrt.(1.0 ./ sum(1 ./ delta_k.^2) )

        mu_limit = (center_mu - k_spread*spread, center_mu + k_spread*spread)
        if !isnothing(mu_lims_dict[pno]) # overwrite mu_lims
            mu_limit = (mu_lims_dict[pno][1], mu_lims_dict[pno][2])
        end 
         
        sig_limit = (max(0, center_sig - k_spread*spread), center_sig + k_spread*spread)
        if !isnothing(sig_lims_dict[pno]) # overwrite sig_lims
            sig_limit = (sig_lims_dict[pno][1], sig_lims_dict[pno][2])
        end 

        # calculate the ranges 
        mu_values = collect(LinRange(mu_limit[1], mu_limit[2], n_points))
        sig_values = collect(LinRange(sig_limit[1], sig_limit[2], n_points))
        p_mu_sig, p_sig, p_mu, n_tot, nn_marg = hyperparamDistTIGER(mu_values, sig_values, dphi0_k, delta_k)
                                                # TODO: right now, using the old version of this distribution
                                                # The new version has some issue I could not fix yet. 
                                                # Maybe someone may try reimplement this function ?
        
        println("Consistency check for calculated distribution:")
        println("Difference in normalization = $(n_tot-nn_marg)")
        
        # run the MCMC
        if MCMC
            println("Running MCMC")
            chain = run_MCMC(dphi0_k, delta_k, mu_limit, sig_limit[2])
            println("MCMC finished")
        end
        
        # create the plot
        println("Injected values: ", val_injected_dict[pno])
        title_str = "PN = $(pno)"
        splot = distributionSummaryPlot(mu_values, sig_values, p_mu_sig, p_sig, p_mu, val_inj=val_injected_dict[pno], title = title_str)
        
        # save plot
        mkpath(figure_dir* nn * "/")
        fig_file_name = figure_dir * nn * "/hyperdist_plot_pn" * pn_order_dic[pno][2] * ".pdf"
        println("\nSaving plot in: $(fig_file_name)")
        savefig(splot, fig_file_name)

        if MCMC
            chain_mu_sigma = [chain[:mu].data, chain[:sigma].data]
            plot_ = plot2DContourKDE(chain_mu_sigma[1], chain_mu_sigma[2], val_injected_dict[pno])
            fig_file_name_MCMC = figure_dir * nn * "/hyperdist_plot_pn_MCMC" * pn_order_dic[pno][2] * ".pdf"
            println("\nSaving MCMC plot in: $(fig_file_name_MCMC)")
            savefig(plot_, fig_file_name_MCMC)

            # save the chain
            serialize( folder * chain_file_name, chain)

            # produce a plot of the contour on top of the distribution

            plot_on_top = plot2DContourKDE!(splot, chain_mu_sigma[1], chain_mu_sigma[2], val_injected_dict[pno])
            fig_file_name_MCMC_on_top = figure_dir * nn * "/hyperdist_plot_pn_MCMC_on_top" * pn_order_dic[pno][2] * ".pdf"
            println("\nSaving MCMC on top plot in: $(fig_file_name_MCMC)")
            savefig(plot_on_top, fig_file_name_MCMC_on_top)

        end

    end
end

