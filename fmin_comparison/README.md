# GW150914-like fmin comparison

This folder reproduces the GW150914-like ET comparison across the low-frequency
cutoffs

```julia
fmin = [2, 5, 10, 15, 20] Hz
```

The workflow ports the old hand-edited `C_conditioned_plot.jl` recipe into an
explicit CairoMakie workflow. The Fisher products are stored as one HDF5 file
per `fmin`, and the plotting script reads those files to build the split
PN-order figure.

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

The default plot overlays the GW150914-like ET `T` injection reference values
used in the old workflow:

```julia
[4.5e-5, 0.0035, 0.011, 0.008, 0.004, 0.034, 0.012, 0.016, 0.088, 0.048]
```

## Configuration notes

The default config uses 100 GW150914-like realizations with fixed
detector-frame chirp mass, symmetric mass ratio, and distance, while spins,
sky position, inclination, polarization, coalescence time, and coalescence
phase are sampled from the broad ranges documented in the TOML file. The plot
uses `events_per_realization = 1`, matching the old fmin-comparison setup.
