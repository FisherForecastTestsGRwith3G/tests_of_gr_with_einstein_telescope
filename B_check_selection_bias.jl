using HDF5
using Plots

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
data_folder_name = user_configs["path_output"]*simulation_tag
output_folder_name =  user_configs["path_output"]*simulation_tag*"/plots/"

# get global indices
global_index_fisher = Dict()
global_index_snr = Dict()
global_index_total = Dict()

for nn in configs["network_list"]

    file_name = data_folder_name *"/data/" * nn *"/global_indices.h5"
    h5open(file_name, "r") do gi_file
        global_index_fisher[nn] = read(gi_file, "fisher")  
        global_index_snr[nn] = read(gi_file, "snr")
        global_index_total[nn] = read(gi_file, "total")
    end

end

# load the catalog with the pn deviations and plot the statistik for potential selection effects
gr_parameter_names = 
    ["eta","chi1","chi2", "dL", "theta", "phi", "iota", "psi", "phiCoal", "tcoal", "lambda1", "lambda2"]

for nn in configs["network_list"]

    file_name = data_folder_name * "/catalog_w_deviations.h5"
    h5open(file_name, "r") do catalog_file
        mkpath(output_folder_name * "/"*nn*"/selection_bias_checks/")        
        param_group = catalog_file["parameter"]
        
        for param_name in keys(param_group)
            
            param_values = read(param_group, param_name)
            param_values_selected = param_values[global_index_total[nn]]
            
            phist = histogram(
                param_values,
                color = :gray,
                alpha = 0.5,
                bins = 30,
                normalized = true,
                label= "Available events"
                )

            phist = histogram!(phist,
                param_values_selected,
                    color = :blue,
                    alpha = 0.5,
                    normalized = true,
                    bins = 30,
                    label= "Selected events")

            xlabel!(phist, param_name)
            #title!(phist, "Difference in distribution for: $(param_name)")
            
            fig_file_name = output_folder_name * "/"*nn*"/selection_bias_checks/"*param_name*".png"
            println("\nSaving plot for parameter '$(param_name)' to: $(fig_file_name)")
            savefig(phist, fig_file_name)
            
        end
    end

    

end