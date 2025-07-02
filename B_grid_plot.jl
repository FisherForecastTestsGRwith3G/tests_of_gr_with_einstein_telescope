using HDF5
using Trapz
using LaTeXStrings
using KernelDensity, Statistics
using Serialization
using Turing
using Base.Threads
using Plots

# include required scripts 
include("_hierarchical_dist.jl")
include("_hierarchical_plotting_utils.jl")
include("_parse_config.jl")
include("_grid_plot_utils.jl")

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 0
    config_file_name = ARGS[1]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForZ(config_file_name)
configs_B, simulation_tag_B = readConfigForB(config_file_name)




mu_vec = configs["mu"]
sigma_vec = configs["sigma"]
pno = configs["pn_waveforms"][1]
PN_name = pn_order_dic[pno][2]
network = configs["network_list"][1]
n_events = configs["n_events"]
n_median = configs["n_median"]
gridSize = length(mu_vec) * length(sigma_vec) 

# set up folder names
data_folder_name = user_configs["path_output"]*simulation_tag * "/"
output_folder_name =  user_configs["path_output"]*simulation_tag*"/results/" * network * "/" *PN_name *"/"
plot_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/" * network * "/" *PN_name *"/"

mkpath(output_folder_name)
mkpath(plot_folder_name)


println("\n Anlayzing grid composed of $gridSize grid points\n")
println("\n Calculating $(n_median) median values for each grid point\n")

# get global index
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()


# process all the plots


idx=0
global res = zeros((length(mu_vec), length(sigma_vec)))
for i in 1:length(mu_vec)
    for j in 1:length(sigma_vec)
        global idx += 1
        println("\n Analyzing data for index $(idx) out of $(gridSize)")
        # skip the (0.,0.) case
        if mu_vec[i] == 0. && sigma_vec[j] == 0.
            continue
        end

        mu = mu_vec[i]
        sigma = sigma_vec[j]
        header = "grid_PN_$(PN_name)_n_$(idx)_mu_" * string(mu) * "_sigma_" * string(sigma)
        data_grid_folder_name = data_folder_name * header * "/data/" * network# * "/pn_" * PN_name 

        println(" Mu: ", mu, " Sigma: ", sigma)
        file_name = data_grid_folder_name * "/global_indices.h5"
        h5open(file_name, "r") do gi_file
            global_index_fisher[network] = read(gi_file, "fisher")  
            global_index_snr[network] = read(gi_file, "snr")
            global_index_total[network] = read(gi_file, "total")
        end

        data_grid_folder_name_pno = data_folder_name * header * "/data/" * network * "/pn_" * PN_name
        filename = data_grid_folder_name_pno * "/single_event_measurements.h5"
        data = Dict()
        h5open(filename, "r") do file
            for key in keys(file)
                data[key] = read(file, key)
            end
        end

        # extract only valid data-points via global index
        dphi0_k = data["dphi0_k"][global_index_total[network]]
        delta_k = data["delta_k"][global_index_total[network]]

        global res[i,j] = reshuffling_bisection(n_median, wrapper_3sigma, dphi0_k, delta_k, mu, sigma)[1]
    end
end

# save data

h5open(output_folder_name*"results.h5", "w") do file
    write(file, "results", res)
end
# else
#     println("The simulation will not be run!")
#     # load the data
#     file_name = output_folder_name * "results.h5"
#     h5open(file_name, "r") do file
#         global res = read(file, "results")
#     end
# end

fig_file_name = plot_folder_name * "grid_plot.pdf"
println("\nSaving grid plot in: $(fig_file_name)")


using CairoMakie

### leave unchanged
fontsize_theme = Theme(fontsize = 24)
set_theme!(fontsize_theme)

MT = Makie.MathTeXEngine
mt_fonts_dir = joinpath(dirname(pathof(MT)), "..", "assets", "fonts", "NewComputerModern")

set_theme!(fonts = (
    regular = joinpath(mt_fonts_dir, "NewCM10-Regular.otf"),
    bold = joinpath(mt_fonts_dir, "NewCM10-Bold.otf")
))
####

tickfontsize = 24
titlefontsize = 26
labelfontsize = 30
lim = collect(0:0.4:3)
lim_label = []
for (i, lim_i) in enumerate(lim)
    if lim_i < 10.

        push!(lim_label, string.(Int.(round.(10 .^lim_i, sigdigits =1))))
    else
        push!(lim_label, string.(Int.(round.(10 .^lim_i, sigdigits =1))))
    end

end

fig = Figure(size = (800, 600))
ax = Axis(fig[1, 1], xlabel = L"\mu", ylabel = L"\sigma", title = "PN: "* string(pno))

# Create a heatmap with the results
hm = Makie.heatmap!(ax, mu_vec, sigma_vec , log10.(res), colormap = :viridis)
Makie.Colorbar(fig[1, 2], hm,  size = 20,
                     ticklabelsize = 20, ticks = (lim, lim_label))

Makie.Label(fig[1, 2, Top()], L"n_{\text{events}}", fontsize = labelfontsize)

ax.titlesize = titlefontsize
ax.xlabelsize = labelfontsize
ax.ylabelsize = labelfontsize

ax.xticklabelsize = tickfontsize
ax.yticklabelsize = tickfontsize
ax.yticklabelspace = 50.


Makie.save(fig_file_name, fig)
