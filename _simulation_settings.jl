include(".paths.jl")

### Set the simulations settings below, to consistently export them in the several other scripts




### Generic settings

PhD = "Matteo"
simulation_tag = "BGR" # how is this simulation called

network_names = ["LHV"]
                #["ETS", "network_0_15km", "network_45_15km", "LHV", "LHVK", "LHV_O3"] # network specs

# specify the PN-orders we want to analyze - may be overriden in the specific scripts
pn_orders = ["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"]   
                #["0", "0.5", "1"] 




### Settings specific to script A_run_catalog_script.jl

needToCreateCatalog = false # set to true if you want to create a new catalog
needToRun = true # set to true if you want to run the simulation
HM = false # set to true if you want to run the simulation with the HM waveform

# specs of catalog
n_events      = 10000
source_type   = "BBH"
catalog_name  = "BGR_TIGER_10k.h5"


# snr threshold
snr_thresh = 12.




### Settings specific to script B_analyze_catalog.jl




### Settings specific to script B_plot_hyperparam_dist.jl

MCMC = true
chain_file_name = "MCMC_chain.jls"
#MCMC_multiple_parallel_chains = true
MCMC_chain_points = 2000
#MCMC_n_chains_multiple_chains = 16




### Additional shared setup

path_catalog, path_output = whoIsThere(PhD)
output_folder_name = path_output*"output/"*simulation_tag*"/"
figure_dir = path_output*"output/"*simulation_tag*"/plots/" # specify where to save the plots

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