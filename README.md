# BGR_with_GWJulia

## Scripts

### Script names

There is two conventions for scripts
* '_nam.jl'  : These scripts are "private" and are called within other scripts. 
* 'L_name.jl : These scripts are to be run from the console. The letters "L" give the order in which to run them.

### Short description

**alpha_generate_catalog_script.jl** <br>
Generates the catalog.  
The number and type of sources to be included in the catalog can be set in the config files.  
Can be run just once, as it is independent of detector networks and subsequent settings.  
The script is called with one argument  

`julia alpha_generate_catalog_script.jl <config_file_name>`

**A_run_catalog_script** <br>
Reads a catalog, creates GR deviations, saves the catalog + GR deviations, evaluates the fisher matrices,
analysis the invertability and snr, stores everything.  
The script is called with two arguments.

`julia A_run_catalog_script.jl <needToEvaluateFisherSNRs> <config_file_name>`

Here: 
* `needToEvaluateFisherSNRs` is either 1 or 0 (by default 1). Can be set to 0 to load precomputed Fisher matrices and SNRs values.
* `config_file_name` is the name of the config file to be used

**B_analyze_catalog.jl** <br>
Visualizes the data in the catalog.  
The script is called with one argument

`julia B_analyze_catalog.jl <config_file_name>`

**B_plot_hyperparam_dist.jl** <br> 
Visualizes the hyperparameter distribution.  
The script is called with one argument.

`julia B_plot_hyperparam_dist.jl <config_file_name>`

**B_check_selection_bias.jl** <br> 
Visualizes biases introduced by the selection criteria due to the invertibility of the fisher.  
The script is called with one argument.

`julia B_check_selection_bias.jl <config_file_name>`

**C_conditioned_plot.jl** <br> 
Produces the plot showing the 90% upper limits for the delta_phi PN deformation coefficients, evaluated from their posterior distributions obtained in hierarchical framework conditioned on sigma = 0.  
From the config files it is possible to:
* show the error bars on the upper limits due to sampling from specific observation realizations (both due to catalogue realization and noise realization).
* show the upper limit from each single events, or the distribution of single events upper limit with violin plots.
* overlay the corresponding LVK GWTC-3 results.

The script is called with one argument.

`julia C_conditioned_plot <config_file_name>`

**D_delta_phi_dist_plot** <br> 
Produces the plot showing the posterior distribution for the delta_phi PN deformation coefficients, obtained in hierarchical framework.  
The plot shows both the generic posterior distribution, and the distribution obtained by conditioning on the hyperparameter sigma = 0.  
The script is called with one argument.

`julia D_delta_phi_dist_plot <config_file_name>`

### How to use config files

The config files are .json files. They specify the parameters of the script and are useful to document different runs.  
For creating such a file see the template file in `/config_files/config_template.json`.

## File structure
The scripts store all their outputs in files structured as follows.

* output/simulation_tag/
    * network/
        * data/network/
            * global_indices.h5 
            * pn_xy/
                * fishers.h5
                * single_event_measurement.h5
                * snrs.h5
                * inspiral_snrs.h5
        * plots/network/
            * pn_xy/
                * hyperdist_plot_pn_MCMC_on_top.pdf
                * hyperdist_plot_pn_MCMC.pdf
                * hyperdist_plot.pdf
                * mudist_plot.pdf
                * catalog_summary/
                    * param_summary_pn_xy.png
                    * snr_error_catalog_pn_xy.png
            * selection_bias_checks/
                * param_name.pdf
    * plots/
        * plot_conditioned_delta_phi_upper_limits_simulation_tag.pdf
        * plot_delta_phi_posterior_dist_simulation_tag.pdf
    * debug/
        * script_D/
            * MCMC_posterior_distribution_pn_xy.pdf

### File-descriptions

**global_indices.h5** <br>
Stores logical arrays, indicating for which events the fisher was invertable for all pn_orders used, which events exceeded the SNR (and eventually inspiral SNR) treshold for all pn_orders, and combinations of such conditions. The key "total" is set to be equivalent to either "total_with_inspiral_snr_cut" or "total_without_inspiral_snr_cut", depending on the config file settings at the time of execution of A_run_catalog_script. <br>
Has the following keys:
* "fisher"
* "snr"
* "inspiral_snr"
* "total" 
* "total_with_inspiral_snr_cut"
* "total_without_inspiral_snr_cut"

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

**inspiral_snrs.h5** <br>
Contains all the inspiral snr values and a logical index for when the events inspiral snr exceeds the set threshold. <br>
Has the following keys:
* "values"
* "index"

**single_event_measurement.h5** <br> 
Has the estiamted error ($\Delta_k$) for the pn-order gr deviation estimated from the fisher, the true parameter-value ($\delta\phi_{Tk}$)  of the gr deviation of the signal and the expected measured value ($\delta\phi_{0k}$). 
Has the following keys:
* "delta_k"
* "dphit_k"
* "dphi0_k" 
