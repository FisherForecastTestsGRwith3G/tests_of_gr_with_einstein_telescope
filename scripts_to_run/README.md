# Tests of GR with HLV

This pipeline generates a BBH catalog, injects beyond-GR post-Newtonian deviations, runs the Fisher analysis for the requested detector networks and waveform families, and stores the resulting summary products in an HDF5 file.

The code can be used to produce Figure 1 of our paper. 

## Reproducing results

Below, we provide instruction on how to rerun scripts in order to produce the plots in our paper. If a script with `A_` or `B_` has already been executed for a plot, it does not need to be run again. 

### Figure 1

```
julia --project=. tests_of_gr_HLV/A_fisher_analysis.jl tests_of_gr_HLV/config_files/config_catalog_200k.toml   
julia --project=. tests_of_gr_HLV/B_population_analysis.jl tests_of_gr_HLV/config_files/config_catalog_200k.toml   
julia --project=. tests_of_gr_HLV/C_fig1.jl tests_of_gr_HLV/config_files/config_catalog_200k.toml   
```

### Figure 2
```
julia --project=. tests_of_gr_HLV/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km0_200k.toml
julia --project=. tests_of_gr_HLV/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml 
julia --project=. tests_of_gr_HLV/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ETS_200k.toml

julia --project=. tests_of_gr_HLV/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ET15km0_200k.toml
julia --project=. tests_of_gr_HLV/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml 
julia --project=. tests_of_gr_HLV/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ETS_200k.

julia --project=. scripts_to_run/C_plot_fig2.jl \
  scripts_to_run/config_files/config_catalog_ET15km0_200k.toml \
  scripts_to_run/config_files/config_catalog_ET15km45_200k.toml \
  scripts_to_run/config_files/config_catalog_ETS_200k.toml
```

### Figure 7

Figure 7 is the CairoMakie rewrite of the old
`old /D_delta_phi_dist_plot.jl` workflow. The old script produced grouped violin
plots for the hierarchical posterior distribution of the PN deformation
coefficients, with an optional outline for the `sigma = 0` conditioned
distribution.

For each configured waveform family and PN order, the rewritten script reads
the Fisher likelihood summaries (`dphi_k`, `delta_k`) for a catalog
realization, or the selected population-analysis summaries when no realization
is requested. It then evaluates the full hierarchical deviation posterior from
Appendix C/Eq. (2.22), drawing that profile directly in CairoMakie. The
conditioned outline is evaluated with `HierDist.getNaiveDistributionOnGrid`,
corresponding to the `sigma = 0` distribution in Eq. (2.26).

```
julia --project=. A_fisher_analysis.jl config_files/config_catalog_ET15km45_200k.toml
julia --project=. B_population_analysis.jl config_files/config_catalog_ET15km45_200k.toml
julia --project=. C_plot_fig7.jl config_files/config_catalog_ET15km45_200k.toml
```

Outputs are written to `results/plots/fig_7/` unless the config overrides
`[plots].outdir`.

For the ET network-comparison version, run the population step for the three
single-network configs first, then run the dedicated plotting config:

```
julia --project=. B_population_analysis.jl config_files/config_catalog_ETS_200k.toml
julia --project=. B_population_analysis.jl config_files/config_catalog_ET15km45_200k.toml
julia --project=. B_population_analysis.jl config_files/config_catalog_ET15km0_200k.toml
julia --project=. C_plot_fig7.jl config_files/config_fig7_ET_200k.toml
```

`config_fig7_ET_200k.toml` overlays the three ET detector networks in the same
panels and uses only `PhenomHM`. Change `[plots.fig7].waveform_family` to
`"PhenomD"` if you want the single-waveform comparison made with IMRPhenomD
instead.

Optional Figure 7 controls can be added under `[plots.fig7]` in the TOML config:

```
[plots.fig7]
grid_points = 700
waveform_family = "PhenomHM"
n_events = 100000
number_of_events_single_realization = 10000
n_events_and_number_of_events_single_realization_refer_directly_to_observed_events = false
realization_index = 1
use_inspiral_snr_threshold = false
plot_conditioned_distribution = true
subplots_pn_order_grouping = [[1], [2, 3, 4, 5], [6, 7, 8, 9, 10]]
offset_x_axis = 0.16
violin_width = 0.34
# y_axis_limits = [[-1e-5, 1e-5], [-0.05, 0.05], [-1.0, 1.0]]
```

### Figure 8

Figure 8 is the CairoMakie rewrite of the old
`old /E_overlayed_hyperparam_dist.jl` plus
`old /python_plots/hyperparameter_contours/B_plot_hyperparam_overlay.py`
workflow. It reads the Fisher-analysis HDF5 output, applies the shared
observation cuts across the requested PN orders, builds one catalog realization,
computes the hierarchical `P(mu, sigma | D)` distribution with
`HierDist.hyperparamDistTIGER`, and overlays the 90% credible contours directly
in Julia.

```
julia --project=. A_fisher_analysis.jl config_files/config_fig8_ET15km45_200k.toml
julia --project=. C_plot_fig8.jl config_files/config_fig8_ET15km45_200k.toml
```

If the ET `2L_45` Fisher file already exists, only the second command is
needed. Outputs are written to `results/plots/fig_8/`.

Optional Figure 8 controls can be added under `[plots.fig8]` in the TOML config:

```
[plots.fig8]
waveform_family = "PhenomHM"
grid_points = 700
credible_interval = 0.90
n_events = 100000
number_of_events_single_realization = 10000
n_events_and_number_of_events_single_realization_refer_directly_to_observed_events = false
realization_index = 1
use_inspiral_snr_threshold = false
minus_one_pn_scale = 1000.0
inset_axis_scale = 10.0
main_x_limits = [-0.010, 0.020]
main_y_limits = [0.0, 0.050]
inset_x_limits = [-0.001, 0.001]
inset_y_limits = [0.0, 0.002]
```

### Figure 9
WIP

## Performing additional checks

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
