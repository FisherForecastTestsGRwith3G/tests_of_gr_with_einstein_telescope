# Grid Plot (Figure 4)

This folder produces the grid plot (Figure 4 in the paper) showing how many observed events are needed
for GR to fall outside the 3-sigma hyperparameter contour at each injected
`(mu, sigma)` point.

The workflow has three stages:

- `A_grid_analysis.jl`: creates one Fisher-analysis config for each grid point
  in `config_files/grid/` and runs the Fisher analysis for those injections.
- `B_grid_plot.jl`: reads the Fisher products, computes the median event count
  needed at each grid point, and writes `results.h5`.
- `C_grid_plot.jl`: reads `results.h5` and saves the final heatmap.

## Reproduce the Plot

Run from the repository root:

```bash
julia --project=. grid_plot/A_grid_analysis.jl grid_plot/config_files/config_grid.toml
julia --project=. grid_plot/B_grid_plot.jl grid_plot/config_files/config_grid.toml
julia --project=. grid_plot/C_grid_plot.jl grid_plot/config_files/config_grid.toml
```

If the per-grid Fisher files already exist, start from `B_grid_plot.jl`. If
`results.h5` already exists, only `C_grid_plot.jl` is needed to redraw the
figure.

## Outputs

The per-grid configs are written to:

```text
grid_plot/config_files/grid/
```

The intermediate grid summary is written to:

```text
grid_plot/results/data/grid/<network>/<pn_order>/results.h5
```

The final figure is written to:

```text
grid_plot/results/plots/grid/<network>/<pn_order>/grid_plot.pdf
```

With the default config, this is:

```text
grid_plot/results/plots/grid/network_45_15km/one/grid_plot.pdf
```

## Configuration

The default configuration is `config_files/config_grid.toml`.

Important fields are:

- `[deviations].pn_orders`: PN order to scan. The plotting scripts expect one
  PN order at a time.
- `[deviations].mu_vec` and `[deviations].sigma_vec`: injected grid values.
- `[catalog].n_events`: number of generated events per grid point.
- `[detectors].network`: detector network name used in the Fisher products.
- `[fisher].waveform`: waveform family. The grid scripts use the first entry.
- `[grid].n_median`: number of reshufflings used to estimate the median event
  count.
