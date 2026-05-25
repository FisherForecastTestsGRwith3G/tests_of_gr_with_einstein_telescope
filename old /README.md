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

`julia C_conditioned_plot.jl <config_file_name>`

This script allows the possibility of overlaying the results from different waveform models and minimum frequency fmin, as explained in the following.

To outline the procedure to overlay different waveform models, we take the plot for the comparison of PhenomD and PhenomHM results in the LVK network as a concrete example. This can be performed by:
* if needed, create a catalog using `julia alpha_generate_catalog_script.jl config_files/config_LVK_final_joined_D_HM_plotC.json`.
* Run `julia A_run_catalog_script.jl 1 config_files/config_LVK_final.json` to evaluate the PhenomHM Fisher results.
* Run `julia A_run_catalog_script.jl 1 config_files/config_LVK_PhenomD_final.json` to evaluate the PhenomD Fisher results. 
* Set `overload_detector_networks_as_waveform_models = true` in the script `C_conditioned_plot.jl`, and appropriately set `list_header_simulation_tags` to contain the value of the variable `header` reported the config files previously used. In particular, for this case, `list_header_simulation_tags = ["final_run_LVK", "final_run_LVK_PhenomD"]`.
* Run `julia C_conditioned_plot.jl config_files/config_LVK_final_joined_D_HM_plotC.json`, answering `y` when prompted.
* Undo the modification implemented in script `C_conditioned_plot.jl` if not needed anymore.

With a similar procedure, it is possible also to create the plot for the GW150914 event, with varying fmin frequency. In this case, the procedure is:
* create first a placeholder catalog structure with (for example) 100 events, running `julia alpha_generate_catalog_script.jl config_files/GW150914_like/config_ET_GW150914_like_fmin_joined_plotC.json`.
* Overload such catalog with several realization of a GW150914-like event. To do so, it is necessary to modify script `A_run_catalog_script.jl`, setting `override_catalog_with_specific_event = true`, and eventually setting the `overriden_` variables to the needed values.
* Run `julia A_run_catalog_script.jl 1 config_files/GW150914_like/_` for all the config files in folder `config_files/GW150914_like` except `config_ET_GW150914_like_fmin_joined_plotC.json`, in order to obtain the values of the Fishers for varying values of fmin. When prompted by the script, answer `y`.
* Set `overload_detector_networks_as_waveform_models = true` in the script `C_conditioned_plot.jl`, and appropriately set `list_header_simulation_tags` to contain the value of the variable `header` reported the config files previously used. In particular, for this case, `list_header_simulation_tags = ["ET_GW150914_like_fmin_2Hz", "ET_GW150914_like_fmin_5Hz","ET_GW150914_like_fmin_10Hz","ET_GW150914_like_fmin_15Hz","ET_GW150914_like_fmin_20Hz"]`.
* Set `override_points_to_be_overlaid = true` in the script `_conditioned_plot_utils.jl`, and appropriately set the variable `LVK_GWTC3_results` to the contain the values of the upper bounds to be shown. For example `LVK_GWTC3_results = [4.5e-5, 0.0035, 0.011, 0.008, 0.004, 0.034, 0.012, 0.016, 0.088, 0.048]` to reproduce the ET bluebook injection.
* Run `julia C_conditioned_plot.jl config_files/GW150914_like/config_ET_GW150914_like_fmin_joined_plotC.json`, answering `y` when prompted.
* Undo the modification implemented in script `A_run_catalog_script.jl`, `C_conditioned_plot.jl` and `_conditioned_plot_utils.jl` if not needed anymore.

When generalizing such procedures to different scenarios, it may be needed to update the `labels_from_networks` function in `_plot_style.jl`.

**C_conditioned_plot_N_dependence.jl** <br> 
Produces the plot showing the trend of the 90% upper limits for the delta_phi PN deformation coefficients (evaluated from their posterior distributions obtained in hierarchical framework conditioned on sigma = 0) as a funcion of the number of (observed) events, for a single chosen detector.  
The config from the config files are mostly overloaded in the script itself: therefore one should read the source code before using it, to make sure the results will be produced as intended.  
The script is called with one argument.

`julia C_conditioned_plot_N_dependence.jl <config_file_name>`

**C_conditioned_plot_single_event_bounds_vs_parameters.jl** <br> 
Produces the plot showing the correlation between the single events 90% upper limits for the delta_phi PN deformation coefficients (evaluated from their posterior distributions obtained in hierarchical framework conditioned on sigma = 0), and other quantities (for example SNR on the y axis, and chirp mass on the z axis, represented with a colorbar).  
The config from the config files are mostly overloaded in the script itself: therefore one should read the source code before using it, to make sure the results will be produced as intended.  
The script is called with one argument.

`julia C_conditioned_plot_single_event_bounds_vs_parameters.jl <config_file_name>`

**D_delta_phi_dist_plot.jl** <br> 
Produces the plot showing the posterior distribution for the delta_phi PN deformation coefficients, obtained in hierarchical framework.  
The plot shows both the generic posterior distribution, and the distribution obtained by conditioning on the hyperparameter sigma = 0.  
The script is called with two arguments.

`julia D_delta_phi_dist_plot.jl <perform_MCMC_sampling> <config_file_name>`

Here: 
* `perform_MCMC_sampling` is either 1 or 0 (by default 1). Can be set to 0 to load precomputed MCMC sampling of the posterior distributions.
* `config_file_name` is the name of the config file to be used

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
        * script_D/network/
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
Has the estiamted error ($\Delta_k$) for the pn-order gr deviation estimated from the fisher, the true parameter-value ($\delta\varphi_{Tk}$)  of the gr deviation of the signal and the expected measured value ($\delta\phi_{0k}$). 
Has the following keys:
* "delta_k"
* "dphit_k"
* "dphi0_k" 
