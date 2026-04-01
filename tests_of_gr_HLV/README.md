# Tests of GR with HLV

This pipeline generates a BBH catalog, injects beyond-GR post-Newtonian deviations, runs the Fisher analysis for the requested detector networks and waveform families, and stores the resulting summary products in an HDF5 file.

The code can be used to produce Figure 1 of our paper. 

## Reproducing results

Below, we provide instruction on how to rerun scripts in order to produce the plots in our paper. If a script with `A_` or `B_` has already been executed for a plot, it does not need to be run again. 

### Figure 1

```
julia --project=. A_fisher_analysis config_files/config_catalog_200k.toml   
julia --project=. B_population_analysis config_files/config_catalog_200k.toml   
julia --project=. C_fig1 config_files/config_catalog_200k.toml   
```

### Figure 9
WIP

## Performin additional checks

### Gaussianity of bootstrap samples in log-space

Checks the distribution of the samples of the 90% upper bounds, calculated from the bootstrapping analysis, resembles a log-normal distribution. Plots the distribution and  

```
julia --project=. A_fisher_analysis config_files/config_catalog_200k.toml   
julia --project=. B_population_analysis config_files/config_catalog_200k.toml   
julia --project=. D_plot_log_normal_check.jl config_files/config_catalog_200k.toml   
```

## Other useful information

### Storage of results and script outputs.
The output folder is always documented in the config files. 
The results of the [fisher analysis](./A_fisher_analysis.jl)) are stored as `.h5` files and their keys to access the data are documented in the docstring of the functions [`write_catalog_to_hdf5()`](./A_fisher_analysis.jl#L31) and [`write_pn_results_to_hdf5()`](./A_fisher_analysis.jl#L74). 
The results of the [population analysis](./B_population_analysis.jl) are stored as `.h5` files and their keys to access the data are documented in the docstring of [`write_population_results_hdf5()`](./B_population_analysis.jl#L76). d
