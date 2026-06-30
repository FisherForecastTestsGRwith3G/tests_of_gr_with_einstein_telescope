# Additional Plots

This directory contains standalone plotting scripts for figures that do not
depend on the main HDF5 analysis pipeline.

## Figure 5: PSD/ASD Curves Used in the Analysis

Figure 5 is produced by [`plot_psd_used.jl`](./plot_psd_used.jl). Run the script
from the repository root:

```bash
julia --project=. additional_plots/plot_psd_used.jl
```

The script reads the detector sensitivity curves from:

- [`create_single_event_datasets/psd_data/hlv_curves/O3b/`](../create_single_event_datasets/psd_data/hlv_curves/O3b/)
- [`create_single_event_datasets/psd_data/et_curves/`](../create_single_event_datasets/psd_data/et_curves/)

It writes the figure to:

```text
additional_plots/results/psds_used.pdf
```

No catalog generation, Fisher analysis, or population-analysis output is
required for this plot.
