# Code Tests of GR with ET

This repository contains Julia workflows for producing the plots used in the
tests-of-GR analyses with HLV and Einstein Telescope detector networks. Run
commands from the repository root with:

```bash
julia --project=.
```

## Plot Workflows

- `scripts_to_run/`: main catalog, population, and plotting scripts for the
  paper figures. See `scripts_to_run/README.md` for the figure-by-figure
  commands.
- `fmin_comparison/`: GW150914-like comparison across low-frequency cutoffs.
  See `fmin_comparison/README.md`.
- `grid_plot/`: injected `(mu, sigma)` grid scan and heatmap of the number of
  events needed to exclude GR at 3 sigma. See `grid_plot/README.md`.
- `additional_plots/`: standalone helper plots, currently including detector
  ASD comparisons.
- `create_single_event_datasets/`: catalog and detector-setup utilities used by
  the analysis scripts.
- `hierachical_combination/`: shared hierarchical-distribution and plotting
  utilities used by the plot workflows.


