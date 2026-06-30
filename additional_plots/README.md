# Additional Plots

This directory contains the helper script for reproducing the HLV amplitude
spectral density comparison plot.

## HLV ASD Comparison

Run the script from the repository root:

```bash
cd /home/joachim-pomper/Desktop/dumpyard/tests_of_gr_with_einstein_telescope
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. additional_plots/plot_hlv_asd_comparisons.jl
```

The script reads the H1, L1, and V1 ASD curves for O3a, O3b, and the O4
pre-observing estimates from:

```text
create_single_event_datasets/psd_data/hlv_curves/
```

It writes both PNG and PDF outputs to:

```text
additional_plots/results/hlv_asd_comparisons.png
additional_plots/results/hlv_asd_comparisons.pdf
```

The output directory is created automatically if it does not already exist.
