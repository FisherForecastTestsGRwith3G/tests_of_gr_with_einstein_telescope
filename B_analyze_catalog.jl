using HDF5
using Plots
using Contour
using Trapz

################################################################################
## Specify simulation specs in this part of the script
# ||||||||
# vvvvvvvv

include("_simulation_settings.jl")

# Below here you may ovveride the simulation settings specified in the _simulation_settings.jl file
# Yet for consistency across different scripts it is recommended to set the simulation settings in the _simulation_settings.jl file

# e.g.:
# specify the pn-orders we want to analze
#pn_orders =  ["-1", "0", "0.5", "1"] 
           # ["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"]  


### Restrict_catalog 
n_events = 10000  # number of events from the catalog used     

output_folder_name = output_folder_name*"data/"

# ^^^^^^^^
# ||||||||
## Specify simulation specs in this part of the script
################################################################################


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

### analyze catalog for each selected pn order
event_ks = collect(1:n_events)
for nn in network_names
    
    idx_g_total = global_index_total[nn][1:n_events]
    idx_g_fisher = global_index_fisher[nn][1:n_events]
    mkpath(figure_dir* nn * "/catalog_summary/")

    for pno in pn_orders
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
        psnr = scatter(
            event_ks,
            snr_values,
            label = "event available", 
            mc=:red, 
            ms=2, 
            ma=0.5)
        plot!(psnr, [1, n_events], [snr_thresh,snr_thresh], label = "")
        scatter!(psnr,
            event_ks[idx_g_total],
            snr_values[idx_g_total],
            label = "event used", 
            markershape=:star5,
            mc=:blue, 
            ms=4, 
            ma=0.5
            )
        ylabel!(psnr, "SNR")
        xlabel!(psnr, "k")
        
        # plot fisher error summary
        fisher_errors = data["delta_k"][1:n_events]
        pfis = scatter(
            event_ks,
            fisher_errors,
            mc=:red, 
            ms=2, 
            ma=0.5,
            label = "")
        scatter!(pfis,
            event_ks[idx_g_total],
            fisher_errors[idx_g_total],
            markershape=:star5,
            mc=:blue, 
            ms=4, 
            ma=0.5,
            label = "",
            ticks = :native, 
            xtickfontcolor = RGBA(0,0,0,0))
        ylabel!(pfis, "Δ_k")
        xlabel!(pfis, "k")

        # plot parameter summery 
        l = @layout [ jeff{1.0w, 0.5h} ; hubert{1.0w, 0.5h}]
     
        snr_fish_plot = plot(pfis, psnr, layout = l)
        fig_file_name = figure_dir * nn * "/catalog_summary/snr_error_catalog_" * pno_name * ".png"
        println("\nSaving plot in: $(fig_file_name)")
        savefig(snr_fish_plot, fig_file_name)

        # summary plot of the rest of the data

        dphi0_k = data["dphi0_k"][1:n_events]
        #dphi0_k = data["dphi0_k"][1:n_events]
        psum = scatter(
            dphi0_k,
            fisher_errors,
            mc=:red, 
            ms=2, 
            ma=0.5,
            label = "events")
        scatter!(psum,
            dphi0_k[idx_g_total],
            fisher_errors[idx_g_total],
            markershape=:star5,
            mc=:blue, 
            ms=4, 
            ma=0.5,
            label = "events used"
            #ticks = :native, 
            #xtickfontcolor = RGBA(0,0,0,0)
        )
        ylabel!(psum, "Δ_k")
        xlabel!(psum, "ɸ_0k")
        
        phist_1 = histogram(
            fisher_errors,
            color = :red,
            alpha = 0.5,
            normalized = true,
            label= "")
        histogram!(phist_1,
            fisher_errors[idx_g_total], 
            color = :blue,
            alpha = 0.5,
            normalized = true, 
            label= "")
        xlabel!(phist_1, "Δ_k")
        
        phist_2 = histogram(
            dphi0_k,
            color = :red,
            alpha = 0.5,
            normalized = true,
            label= ""
        )
        histogram!(phist_2,
            dphi0_k[idx_g_total], 
            color = :blue,
            alpha = 0.5,
            normalized = true, 
            label= "")
        xlabel!(phist_2, "ɸ_0k")

        l = @layout [ jeff{0.6w, 0.6h} heisenberg{0.4w, 0.6h} ; dirac{0.6w, 0.4h} david{0.4w, 0.4h}]
        datasum_plot = plot(psum, phist_1, phist_2, layout = l)
        fig_file_name = figure_dir * nn * "/catalog_summary/param_summary_" * pno_name * ".png"
        println("\nSaving plot in: $(fig_file_name)")
        savefig(datasum_plot, fig_file_name)
    end 
end