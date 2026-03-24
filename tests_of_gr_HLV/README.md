# GWjulia BGR Pipeline

This pipeline generates a BBH catalog, injects beyond-GR post-Newtonian deviations, runs the Fisher analysis for the requested detector networks and waveform families, and stores the resulting summary products in an HDF5 file.

## Storage structure

```text
fisher_results_<tag>.h5
├── bbh_catalog
│   ├── .attrs
│   └── parameter
└── <pn_order_name>
    ├── .attrs
    ├── delta_phi_pn
    └── <network>
        └── <waveform_family>
            ├── snr
            ├── isnr
            ├── invc
            ├── delta_k
            └── dphi_k
```

Example dataset path:
`/minus_one/LHV/PhenomHM/snr`

## Config expectations

The analysis expects a TOML config file with the following sections:

- `[general]`
- `[catalog]`
- `[deviations]`
- `[detectors]`
- `[fisher]`

Important fields:

- `catalog.n_events`: Number of events required in the final catalog.
- `catalog.z_cut`: Optional redshift cut. If omitted, it defaults to `Inf`.
- `deviations.pn_orders`: List of PN orders to analyze.
- `deviations.mu`: List of Gaussian means, one per PN order.
- `deviations.sigma`: List of Gaussian standard deviations, one per PN order.
- `detectors.network`: List of detector networks.
- `fisher.waveform`: Must be given as a list, for example `["PhenomHM"]`.
- `fisher.fmin`: Lower frequency cutoff used in the Fisher analysis.

## Attributes

`/bbh_catalog/.attrs`

- `catalog_tag`: Tag of the simulated catalog used in the run.
- `requested_n_events`: Number of events requested in the config file.
- `seed_initial`: First random seed used to generate the catalog.
- `seed_final`: Last random seed used when extending the catalog after the redshift cut.
- `z_cut`: Redshift threshold applied to the catalog. If no cut is requested, this is `Inf`.

`/<pn_order_name>/.attrs`

- `mu`: Mean of the Gaussian population distribution used to draw the beyond-GR deviation for this PN order.
- `sigma`: Standard deviation of the Gaussian population distribution used to draw the beyond-GR deviation for this PN order.

## Stored datasets

The HDF5 file stores the full BBH catalog used in the Fisher study together with the injected beyond-GR deviations and the derived Fisher-analysis summary quantities for every requested PN order, detector network, and waveform family.

### Catalog data

`/bbh_catalog/parameter` contains the source catalog parameters for the events that actually enter the analysis after applying the redshift cut and any required catalog extension.

### PN-order data

`/<pn_order_name>/delta_phi_pn` contains the injected beyond-GR deviation values for that PN order, one value per event in the final catalog.

### Per network / waveform outputs

`/<pn_order_name>/<network>/<waveform_family>/snr` stores the full matched-filter signal-to-noise ratio for each event.

`/<pn_order_name>/<network>/<waveform_family>/isnr` stores the inspiral-only signal-to-noise ratio for each event.

`/<pn_order_name>/<network>/<waveform_family>/invc` is a Boolean array indicating whether the Fisher matrix inversion was successful for each event.

`/<pn_order_name>/<network>/<waveform_family>/delta_k` stores the Fisher-estimated 1-sigma uncertainty on the deviation parameter for each event.

`/<pn_order_name>/<network>/<waveform_family>/dphi_k` stores the expected measured deviation value for each event after combining the injected deviation with the Fisher uncertainty model.
