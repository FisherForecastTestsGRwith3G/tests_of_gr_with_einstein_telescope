# BGR_with_GWJulia

## Scripts

### Script names

There is two conventions for scripts
* '_nam.jl'  : These scripts are "private" and are called within other scripts. 
* 'L_name.jl : These scripts are to be run from the console. The letters "L" give the order in which to run them.

### Short description

**A_run_catalog_script** <br>
Reads a catalog, creates GR deviations, saves the catalog + GR deviations, evaluates the fisher matrices,
analysis the invertability and snr, stores everything.

**B_analyze_catalog.jl** <br>
Visualizes the data in the catalog.

**B_plot_hyperparam_dist.jl** <br> 
Visualizes the hyperparameter distribution

## File strucutre
The scripts store all their outputs in files structured as follows.

* output/simulation_tag/network/
    * data/network/
        * global_indicies.h5 
        * pn_xy/
            * fishers.h5
            * single_event_measurement.h5
            * snrs.h5
    * plots/network/
        * catalog_summary/
            * param_summary_pn_xy.png
            * snr_error_catalog_pn_xy.png
        * hyperdist_plot_pn_xy.png 

### File-descriptions

**global_indicies.h5** <br>
Stores logical arrays, indicating for which events the fisher was invertable for all pn_orders used, which events exceeded the treshhold for all pn_orders or both. <br>
Has the following keys:
* "fisher"
* "snr"
* "total" 

**fishers.h5** <br>
Contains the fisher matrices and a logical index for the indices whos fisher is invertible. <br>
Has the following keys:
* "values"
* "index"

**snrs.h5** <br>
Contains all the snr values and a logical index for when the events snr exceeds the set threshold. <br>
Has the following keys:
* "values"
* "index"

**single_event_measurement.h5** <br> 
Has the estiamted error ($\Delta_k$) for the pn-order gr deviation estimated from the fisher, the true parameter-value ($\delta\phi_{Tk}$)  of the gr deviation of the signal and the expected measured value ($\delta\phi_{0k}$). 
Has the following keys:
* "delta_k"
* "dphit_k"
* "dphi0_k" 
