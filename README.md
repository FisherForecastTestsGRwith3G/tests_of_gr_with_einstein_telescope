# Code Tests of GR with ET

This repository contains Julia workflows for producing the results of
https://arxiv.org/abs/2511.07520 and plots for HLV and Einstein Telescope
detector networks.

## Citing Our Work

If you use this code or the associated results, please cite:

```bibtex
@article{Begnoni:2025mtz,
    author = "Begnoni, Andrea and Del Pozzo, Walter and Pegorin, Matteo and Pomper, Joachim and Ricciardone, Angelo",
    title = "{Tests of General Relativity with Einstein Telescope}",
    eprint = "2511.07520",
    archivePrefix = "arXiv",
    primaryClass = "gr-qc",
    month = "11",
    year = "2025"
}
```

The INSPIRE record is available at https://inspirehep.net/literature/3081766.

See [AUTHORS.md](AUTHORS.md) for the author list. This repository is released
under the MIT license; see [LICENSE](LICENSE).

## Workflow

Run every command below from the repository root:

```bash
cd path/to/this/repo/tests_of_gr_with_einstein_telescope
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The production configurations are intentionally large. The main catalog configs
generate 200,000 events and some population/scaling steps bootstrap 10,000
catalogs, so expect the `A_` and `B_` stages to be the expensive parts. Plotting
scripts are usually quick once their HDF5 inputs exist.

## Repository Layout

- `scripts_to_run/`: main catalog, Fisher, population, scaling, and plotting
  scripts for Figures 1, 2, 3, 7, 8, and 9.
- `fmin_comparison/`: GW150914-like ET comparison across low-frequency cutoffs.
- `grid_plot/`: injected `(mu, sigma)` grid scan and heatmap of the number of
  events needed to exclude GR at 3 sigma.
- `additional_plots/`: standalone helper plots, currently the HLV ASD
  comparison.
- `create_single_event_datasets/`: catalog and detector setup utilities.
- `hierachical_combination/`: hierarchical-distribution and likelihood
  utilities used by the plot workflows.

## Analysis Stages

The main scripts use a resumable naming convention:

- `A_fisher_analysis.jl` creates catalog and Fisher products:
  `scripts_to_run/results/data/fisher_results_<catalog_tag>.h5`.
- `B_population_analysis.jl` reads Fisher products and writes population
  summaries:
  `scripts_to_run/results/data/population_results_<bootstrap_tag>.h5`.
- `B_scaling_analysis.jl` reads Fisher products and writes scaling summaries:
  `scripts_to_run/results/data/scaling_results_<scaling_tag>.h5`.
- `C_plot_*.jl` reads existing HDF5 products and writes `.png` and `.pdf`
  figures below `scripts_to_run/results/plots/`.
- `D_*.jl` scripts are optional diagnostics and distribution checks.

If a required `.h5` file already exists, you can skip the earlier stage and
rerun only the downstream analysis or plotting script.

## Main Paper Plots

### Figure 1: HLV O3b Single-Event and Population Constraints

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
julia --project=. scripts_to_run/B_population_analysis.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
julia --project=. scripts_to_run/C_plot_fig1.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
```

Output:

```text
scripts_to_run/results/plots/fig_1/fig1_hlv_o3b_bbh_200k_stb_9.{png,pdf}
```

### Figure 2: ET Network Comparison

Run the Fisher and population stages for all three ET network configs:

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km0_200k.toml
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ETS_200k.toml

julia --project=. scripts_to_run/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ET15km0_200k.toml
julia --project=. scripts_to_run/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
julia --project=. scripts_to_run/B_population_analysis.jl scripts_to_run/config_files/config_catalog_ETS_200k.toml

julia --project=. scripts_to_run/C_plot_fig2.jl \
  scripts_to_run/config_files/config_catalog_ET15km0_200k.toml \
  scripts_to_run/config_files/config_catalog_ET15km45_200k.toml \
  scripts_to_run/config_files/config_catalog_ETS_200k.toml
```

Output:

```text
scripts_to_run/results/plots/fig_2/fig2_detector_comparison_ET15km0_bbh_200k_bts_10k.{png,pdf}
```

### Figure 3: Constraint Scaling with Number of Observed Events

The default Figure 3 workflow uses the ET `2L_45` config and the scaling stage:

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
julia --project=. scripts_to_run/B_scaling_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
julia --project=. scripts_to_run/C_plot_fig3.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
```

Output:

```text
scripts_to_run/results/plots/fig_3/fig3_ET15km45_bbh_200k_30_10k.{png,pdf}
```

### Figure 7: Hierarchical PN-Deviation Posteriors

For the HLV O3b version, use the HLV population product:

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
julia --project=. scripts_to_run/B_population_analysis.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
julia --project=. scripts_to_run/C_plot_fig7.jl scripts_to_run/config_files/config_fig7_hlv_o3b_200k.toml
```

For the ET network-comparison version, the config draws fixed-size realizations
directly from the three Fisher products:

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ETS_200k.toml
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_ET15km0_200k.toml
julia --project=. scripts_to_run/C_plot_fig7.jl scripts_to_run/config_files/config_fig7_ET_200k.toml
```

Outputs:

```text
scripts_to_run/results/plots/fig_7/fig7_fig7_hlv_o3b_bbh_200k_stb_9.{png,pdf}
scripts_to_run/results/plots/fig_7/fig7_fig7_ET_bbh_200k_bts_10k_PhenomHM.{png,pdf}
```

### Figure 8: Hyperparameter Contours for One ET Realization

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_fig8_ET15km45_200k.toml
julia --project=. scripts_to_run/C_plot_fig8.jl scripts_to_run/config_files/config_fig8_ET15km45_200k.toml
```

If `scripts_to_run/results/data/fisher_results_ET15km45_bbh_200k.h5` already
exists, only the plotting command is needed.

Output:

```text
scripts_to_run/results/plots/fig_8/fig8_ET15km45_bbh_200k_PhenomHM_realization_1.{png,pdf}
```

### Figure 9: HLV Observability and Fisher-Improvement Diagnostics

Figure 9 reads the HLV Fisher product directly:

```bash
julia --project=. scripts_to_run/A_fisher_analysis.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
julia --project=. scripts_to_run/C_plot_fig9.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
```

Output:

```text
scripts_to_run/results/plots/fig_9/fig9_hlv_o3b_bbh_200k_stb_9.{png,pdf}
```

## Additional Plot Workflows

### GW150914-Like fmin Comparison

This workflow compares ET constraints for:

```julia
fmin = [2, 5, 10, 15, 20] Hz
```

Run:

```bash
julia --project=. fmin_comparison/A_fisher_analysis.jl fmin_comparison/config_files/config_gw150914_fmin_comparison.toml
julia --project=. fmin_comparison/C_plot_fmin_comparison.jl fmin_comparison/config_files/config_gw150914_fmin_comparison.toml
```

Outputs:

```text
fmin_comparison/results/data/fisher_results_gw150914_like_fmin_<value>.h5
fmin_comparison/results/plots/fmin_comparison/fmin_comparison_gw150914_like.{png,pdf}
```

### Injected `(mu, sigma)` Grid Plot

The grid workflow has three stages:

- `grid_plot/A_grid_analysis.jl`: creates one generated Fisher config for each
  grid point and runs the Fisher analysis.
- `grid_plot/B_grid_plot.jl`: reads the per-grid Fisher products and writes
  `results.h5`.
- `grid_plot/C_grid_plot.jl`: reads `results.h5` and saves the heatmap.

For the `+1` PN-order grid:

```bash
julia --project=. grid_plot/A_grid_analysis.jl grid_plot/config_files/config_grid_one.toml
julia --project=. grid_plot/B_grid_plot.jl grid_plot/config_files/config_grid_one.toml
julia --project=. grid_plot/C_grid_plot.jl grid_plot/config_files/config_grid_one.toml
```

For the `-1` PN-order grid:

```bash
julia --project=. grid_plot/A_grid_analysis.jl grid_plot/config_files/config_grid_minus_one.toml
julia --project=. grid_plot/B_grid_plot.jl grid_plot/config_files/config_grid_minus_one.toml
julia --project=. grid_plot/C_grid_plot.jl grid_plot/config_files/config_grid_minus_one.toml
```

Outputs:

```text
grid_plot/config_files/grid/<network>_<waveform>/
grid_plot/results/data/grid/<network>/<pn_order>/results.h5
grid_plot/results/plots/grid/<network>/<pn_order>/grid_plot.pdf
```

If the per-grid Fisher files already exist, start from `B_grid_plot.jl`. If
`results.h5` already exists, run only `C_grid_plot.jl`.

### HLV ASD Comparison

```bash
julia --project=. additional_plots/plot_hlv_asd_comparisons.jl
```

Output:

```text
additional_plots/results/hlv_asd_comparisons.{png,pdf}
```

## Optional Diagnostics

### Bootstrap Log-Normal Check

After running `A_fisher_analysis.jl` and `B_population_analysis.jl` for a config:

```bash
julia --project=. scripts_to_run/D_plot_log_normal_check.jl scripts_to_run/config_files/config_catalog_hlv_o3b_200k.toml
```

Outputs are written to:

```text
scripts_to_run/results/plots/checks/bootstrap/
```

### Hyperparameter Distribution Check

This diagnostic reads a Fisher product, selects `N` events, and evaluates
per-PN-order hyperparameter distributions for one waveform family:

```bash
julia --project=. scripts_to_run/D_hyperparmDist.jl scripts_to_run/config_files/config_catalog_ET15km45_200k.toml 10000 PhenomHM
```

The optional second and third arguments are `N` and `waveform_family`. If they
are omitted, the script uses the config's `bootstrap.n_catalog` and the first
configured waveform family.

Outputs are written below:

```text
scripts_to_run/results/plots/distributions/
```

## Editing Configurations

The TOML files in `scripts_to_run/config_files/`,
`fmin_comparison/config_files/`, and `grid_plot/config_files/` control:

- catalog size and random seed in `[catalog]` and `[general]`;
- PN orders and injected deviations in `[deviations]`;
- detector network in `[detectors]`;
- waveform family and low-frequency cutoff in `[fisher]`;
- SNR cuts in `[population]`;
- bootstrap and scaling sizes in `[bootstrap]` and `[scaling]`;
- output tags and plot-specific settings in `[plots]`.

Most output paths are relative to the directory containing the script. For
example, `outdir = "./results/data/"` in a `scripts_to_run` config writes under
`scripts_to_run/results/data/`.
