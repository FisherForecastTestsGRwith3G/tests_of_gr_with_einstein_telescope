using Pkg
using HDF5
import JSON

parentdir = dirname(pwd())
Pkg.activate(string(parentdir)*"/GW.jl")
using GW

# include required scripts 
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

configs, simulation_tag = readConfigForalpha(config_file_name)

println("Creating catalog with ", configs["n_events"], " number of events and source type: ", configs["source_type"])

# Generate a catalog
@time GenerateCatalog(
    configs["n_events"], 
    configs["source_type"], 
    name_catalog=configs["catalog_name"],
    folder=user_configs["path_catalog"]
)