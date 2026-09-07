# GW150914-like fmin comparison

This folder reproduces the GW150914-like ET comparison across the low-frequency
cutoffs

```julia
fmin = [2, 5, 10, 15, 20] Hz
```

## Scripts

- `A_fisher_analysis.jl`: builds a GW150914-like event ensemble and evaluates
  the Fisher products for each configured `fmin`.
- `C_plot_fmin_comparison.jl`: reads the Fisher products, applies the common
  SNR/Fisher selection across PN orders, computes the conditioned
  `sigma = 0` upper limits, and saves the comparison figure.
- `_config_parser.jl`: shared TOML parser and `fmin` tag helpers.

## Reproduce the figure

Run from the repository root:

```bash
julia --project=. fmin_comparison/A_fisher_analysis.jl fmin_comparison/config_files/config_gw150914_fmin_comparison.toml
julia --project=. fmin_comparison/C_plot_fmin_comparison.jl fmin_comparison/config_files/config_gw150914_fmin_comparison.toml
```

The Fisher files are written to:

```text
fmin_comparison/results/data/
```

The figure is written to:

```text
fmin_comparison/results/plots/fmin_comparison/
```

## Configuration notes

The default config uses 100 GW150914-like realizations with fixed
detector-frame chirp mass, symmetric mass ratio, and distance, while spins,
sky position, inclination, polarization, coalescence time, and coalescence
phase are sampled from the broad ranges documented in the TOML file. The plot
uses `events_per_realization = 1`.
