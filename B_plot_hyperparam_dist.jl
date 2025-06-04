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
include("_setup_MCMC.jl")
include("_parse_config.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForB(config_file_name)

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

# process all the plots
 for nn in configs["network_list"]

    for pno in configs["pn_waveforms"]

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
        n_events_used = configs["n_events"]
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

        ### calculate the hyper-parameter distribution
        # first estimate where to place it
        center_mu = sum(dphi0_k) / n_events_used
        center_sig = max(0, sqrt.(sum((center_mu .- dphi0_k).^2) ./  n_events_used))
        spread = sqrt.(1.0 ./ sum(1 ./ delta_k.^2) )

        k_spread = configs["k_spread"]
        mu_limit = (center_mu - k_spread*spread, center_mu + k_spread*spread)
        if !isnothing(configs["mu_lims"][pno]) # overwrite mu_lims
            mu_limit = (configs["mu_lims"][pno][1], configs["mu_lims"][pno][2])
        end 
         
        sig_limit = (max(0, center_sig - k_spread*spread), center_sig + k_spread*spread)
        if !isnothing(configs["sigma_lims"][pno]) # overwrite sig_lims
            sig_limit = (configs["sigma_lims"][pno][1], configs["sigma_lims"][pno][2])
        end 

        # calculate the ranges 
        mu_values = collect(LinRange(mu_limit[1], mu_limit[2], configs["n_points"]))
        sig_values = collect(LinRange(sig_limit[1], sig_limit[2], configs["n_points"]))
        p_mu_sig, p_sig, p_mu, n_tot, nn_marg = hyperparamDistTIGER(mu_values, sig_values, dphi0_k, delta_k)
        
        println("Consistency check for calculated distribution:")
        println("Difference in normalization = $(n_tot-nn_marg)")
        
        # run the MCMC
        if configs["run_mcmc"]
            println("Running MCMC")
            chain = run_MCMC(dphi0_k, delta_k, mu_limit, sig_limit[2], configs["mcmc_samples"])
            println("MCMC finished")
        end
        
        # create the plot
        println("Injected values: ", configs["injected_values"][pno])
        title_str = "PN = $(pno)"
        splot = distributionSummaryPlot(mu_values, sig_values, p_mu_sig, p_sig, p_mu, val_inj=configs["injected_values"][pno], title = title_str)
        
        # save plot
        pno_name = "pn_" *pn_order_dic[pno][2]
        mkpath(output_folder_name* nn * "/"*pno_name)
        fig_file_name = output_folder_name * nn * "/"*pno_name*"/hyperdist_plot.pdf"
        println("\nSaving hyperparameter distribution plot in: $(fig_file_name)")
        savefig(splot, fig_file_name)

        if configs["run_mcmc"]
            chain_mu_sigma = [chain[:mu].data, chain[:sigma].data]
            plot_ = plot2DContourKDE(chain_mu_sigma[1], chain_mu_sigma[2], configs["injected_values"][pno])
            fig_file_name_MCMC = output_folder_name * nn * "/"*pno_name*"/hyperdist_plot_pn_MCMC.pdf"
            println("\nSaving MCMC plot in: $(fig_file_name_MCMC)")
            savefig(plot_, fig_file_name_MCMC)

            # save the chain
            #TODO: Specify better the chain file name
            serialize( folder * "mcmc_chain", chain)

            # produce a plot of the contour on top of the distribution

            plot_on_top = plot2DContourKDE!(splot, chain_mu_sigma[1], chain_mu_sigma[2])
            fig_file_name_MCMC_on_top = output_folder_name * nn * "/"*pno_name*"/hyperdist_plot_pn_MCMC_on_top.pdf"
            println("\nSaving MCMC on top plot in: $(fig_file_name_MCMC)")
            savefig(plot_on_top, fig_file_name_MCMC_on_top)

        end

        ### Comparison plot between naive and full distribution 
        p_mu_naive = naiveMuDistTIGER(mu_values, dphi0_k, delta_k, maximum(p_mu))
        mu_plot = plot(mu_values, p_mu, color=:blue, linewidth = 2, label = "p(μ|D)")
        plot!(mu_plot, mu_values, p_mu_naive, color=:green, linewidth = 2, label = "p(μ| σ=0, D)")
        xlims!(mu_limit)
        xlabel!(mu_plot, "μ")
        ylabel!(mu_plot, "p(μ)")

        fig_file_name = output_folder_name * nn * "/"*pno_name*"/mudist_plot.pdf"
        println("Saving mu-distribution plot in: $(fig_file_name)")
        savefig(mu_plot, fig_file_name)

    end
end

