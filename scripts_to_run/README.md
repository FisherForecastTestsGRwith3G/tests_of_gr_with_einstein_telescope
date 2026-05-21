# Tests of GR with HLV

This pipeline generates a BBH catalog, injects beyond-GR post-Newtonian deviations, runs the Fisher analysis for the requested detector networks and waveform families, and stores the resulting summary products in an HDF5 file.

The code can be used to produce Figure 1 of our paper. 

## Reproducing results

Below, we provide instruction on how to rerun scripts in order to produce the plots in our paper. If a script with `A_` or `B_` has already been executed for a plot, it does not need to be run again. 

### Figure 1

```
julia --project=. tests_of_gr_HLV/A_fisher_analysis tests_of_gr_HLV/config_files/config_catalog_200k.toml   
julia --project=. tests_of_gr_HLV/B_population_analysis tests_of_gr_HLV/config_files/config_catalog_200k.toml   
julia --project=. tests_of_gr_HLV/C_fig1 tests_of_gr_HLV/config_files/config_catalog_200k.toml   
```

### Figure 9
WIP

### Figure 7

Figure 7 is the CairoMakie rewrite of the old
`old /D_delta_phi_dist_plot.jl` workflow. The old script produced grouped violin
plots for the hierarchical posterior distribution of the PN deformation
coefficients, with an optional outline for the `sigma = 0` conditioned
distribution.

The rewritten script uses the population-analysis HDF5 output written by
[`B_population_analysis.jl`](./B_population_analysis.jl). For each configured
waveform family and PN order it reads `selected_dphi_k` and `selected_delta_k`,
builds the hierarchical distribution with `HierDist.hyperparamDistTIGER`,
marginalizes over `sigma` with `HierDist.getDistributionOnGrid`, and draws the
posterior profiles directly in CairoMakie. The conditioned outline is evaluated
with `HierDist.getNaiveDistributionOnGrid`.

```
julia --project=. A_fisher_analysis.jl config_files/config_catalog_ET15km45_200k.toml
julia --project=. B_population_analysis.jl config_files/config_catalog_ET15km45_200k.toml
julia --project=. C_plot_fig7.jl config_files/config_catalog_ET15km45_200k.toml
```

Outputs are written to `results/plots/fig_7/` unless the config overrides
`[plots].outdir`.

Optional Figure 7 controls can be added under `[plots.fig7]` in the TOML config:

```
[plots.fig7]
grid_points = 700
plot_conditioned_distribution = true
subplots_pn_order_grouping = [[1], [2, 3, 4, 5], [6, 7, 8, 9, 10]]
offset_x_axis = 0.16
violin_width = 0.34
# y_axis_limits = [[-1e-5, 1e-5], [-0.05, 0.05], [-1.0, 1.0]]
```

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
