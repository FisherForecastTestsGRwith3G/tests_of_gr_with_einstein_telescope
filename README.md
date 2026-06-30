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

## Reproducing Main Paper Plots

Instructions to reproduce Figures 1,2,3,7,8, and 9 can be found in [`scripts_to_run/README.md`](scripts_to_run/README.md).
Figure 4 can be reproduced following the instructions found in [`grid_plot/README.md`](grid_plot/README.md).
To reproduce Figure 6 look at [`fmin_comparison/README.md`](fmin_comparison/README.md), while instruction for Figure 5 are provided in [`additional_plots/README.md`](additional_plots/README.md).
