using TOML
using HDF5

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#
"""
    write_catalog_to_hdf5(file, catalog, config, initial_seed, final_seed)

Store the generated binary-black-hole catalog in the output HDF5 file.

The function creates the `/bbh_catalog` group, records the catalog-generation
metadata as attributes, and writes the event-wise source parameters into the
`/bbh_catalog/parameter` subgroup. The stored attributes are intended to make
the resulting file self-describing and include:

- `catalog_tag`: label used to identify the catalog in downstream analyses
- `requested_n_events`: target number of retained events after cuts
- `seed_initial`: first random seed used to build the catalog
- `seed_final`: last random seed needed to reach the requested sample size
- `z_cut`: redshift threshold applied before truncation

The parameter subgroup contains one dataset per catalog field and preserves the
event ordering used throughout the Fisher-analysis pipeline. The dataset keys
written under `/bbh_catalog/parameter` are `mc`, `eta`, `chi_1`, `chi_2`, `dL`,
`theta`, `phi`, `iota`, `psi`, `t_coal`, `phi_coal`, and `z`.
"""
function write_catalog_to_hdf5(file, catalog, config, initial_seed, final_seed)
    catalog_grp = create_group(file, "bbh_catalog")
    attrs(catalog_grp)["catalog_tag"]        = config["catalog_tag"]
    attrs(catalog_grp)["requested_n_events"] = config["n_events"]
    attrs(catalog_grp)["seed_initial"]       = initial_seed
    attrs(catalog_grp)["seed_final"]         = final_seed
    attrs(catalog_grp)["z_cut"]              = config["z_cut"]

    param_grp = create_group(catalog_grp, "parameter")
    write(param_grp, "mc"      , catalog.mc)
    write(param_grp, "eta"     , catalog.eta)
    write(param_grp, "chi_1"   , catalog.chi_1)
    write(param_grp, "chi_2"   , catalog.chi_2)
    write(param_grp, "dL"      , catalog.dL)
    write(param_grp, "theta"   , catalog.theta)
    write(param_grp, "phi"     , catalog.phi)
    write(param_grp, "iota"    , catalog.iota)
    write(param_grp, "psi"     , catalog.psi)
    write(param_grp, "t_coal"  , catalog.t_coal)
    write(param_grp, "phi_coal", catalog.phi_coal)
    write(param_grp, "z"       , catalog.z)
end

"""
    write_pn_results_to_hdf5(file, pno, mu, sigma, pn_deviation, network, wf_fam,
                             snr, isnr, invc, delta_k, dphi_k)

Write the Fisher-analysis products associated with one post-Newtonian deviation
order and one waveform family.

Results are organized hierarchically as:

- `/<pn_order>`: stores the deformation hyper-parameters `mu` and `sigma`
- `/<pn_order>/delta_phi_pn`: stores the injected PN deviation for each event
- `/<pn_order>/<network>/<waveform>`: stores the event-wise analysis outputs

The waveform-level datasets contain the full-network signal-to-noise ratio
(`snr`), inspiral-only signal-to-noise ratio (`isnr`), the Fisher-matrix
invertibility flag (`invc`), the estimated 1-sigma uncertainty on the PN
deformation parameter (`delta_k`), and the injected phase coefficient
(`dphi_k`). Existing PN-order or network groups are reused so that multiple
waveform families can be appended without rewriting shared data.
"""
function write_pn_results_to_hdf5(file, pno, mu, sigma, pn_deviation, network, 
    wf_fam, snr, isnr, invc, delta_k, dphi_k)
    
    pnos = createSED.pnoString(pno)
    pno_grp = haskey(file, pnos) ? file[pnos] : create_group(file, pnos)
    attrs(pno_grp)["mu"]    = mu
    attrs(pno_grp)["sigma"] = sigma
    if !haskey(pno_grp, "delta_phi_pn")
        write(pno_grp , "delta_phi_pn", pn_deviation)
    end

    network_grp = haskey(pno_grp, network) ? pno_grp[network] : create_group(pno_grp, network)
    wf_grp = haskey(network_grp, wf_fam) ? network_grp[wf_fam] : create_group(network_grp, wf_fam)

    write(wf_grp, "snr", snr)
    write(wf_grp, "isnr", isnr)
    write(wf_grp, "invc", invc)
    write(wf_grp, "delta_k", delta_k)
    write(wf_grp, "dphi_k", dphi_k)
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
"""
    run_fisher_analysis(config::Dict)

Execute the full Fisher-analysis stage for the HLV tests-of-GR workflow.

The routine performs four steps:

1. Build a binary-black-hole catalog from the configured random seed.
2. Apply the redshift cut and extend the catalog with additional seeds until
   the requested number of surviving events is available.
3. For each requested waveform family, compute the baseline SNR and inspiral
   SNR using the reference PN order.
4. For every configured PN order, generate synthetic deviations, evaluate the
   single-event Fisher analysis, and write all products to an HDF5 file.

Expected configuration entries include the catalog size and tag, redshift cut,
random seed, PN orders, deviation hyper-parameters (`mu`, `sigma`), detector
network name, waveform families, lower frequency cutoff, SNR thresholds, and
the output directory. The function returns the absolute path of the generated
HDF5 file.

The output file contains both the catalog that was analyzed and all Fisher
products needed by later population-level and plotting scripts, allowing the
downstream stages to run without regenerating the single-event analysis.
"""
function run_fisher_analysis(config::Dict)
    println("Initializing BBH catalog")

    target_n_events = config["n_events"]
    initial_seed = config["seed"]
    next_seed = initial_seed
    output_file = joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")

    #--------------------------------------------------------------------------#
    # Create catalog
    #--------------------------------------------------------------------------#
    catalog = createSED.BBHCatalog(target_n_events, next_seed)
    catalog = createSED.applyRedshiftCut(catalog, config["z_cut"])

    while createSED.length(catalog) < target_n_events
        next_seed += 1
        println(
            "Catalog has $(createSED.length(catalog)) events after z_cut=$(config["z_cut"]). " *
            "Extending with seed $(next_seed)."
        )

        additional_catalog = createSED.BBHCatalog(target_n_events, next_seed)
        additional_catalog = createSED.applyRedshiftCut(additional_catalog, config["z_cut"])
        catalog = createSED.merge(catalog, additional_catalog)
    end

    catalog = createSED.truncateCatalog(catalog, target_n_events)
    println("Catalog initialized with $(createSED.length(catalog)) events.")
    mkpath(dirname(output_file))

    #--------------------------------------------------------------------------#
    # Perform fisher analysis for various PN-orders, networks and waveforms    #
    #--------------------------------------------------------------------------#
    println("Writing results to $(output_file)")
    h5open(output_file, "w") do file
        write_catalog_to_hdf5(file, catalog, config, initial_seed, next_seed)

        reference_pno = first(config["pn_orders"])

        for wf_fam in config["waveform_families"]
            println("\nProcessing waveform family: $(wf_fam)")
            println("Preparing SNR and inspiral SNR calculations using reference PN order $(reference_pno)")

            # Calculate SNR and inspiral SNR first, as we are always evaluating 
            # the waveform in the point where the deviation vanishes.
            snr, isnr = createSED.computeSNRsFromCatalog(
                catalog          ,
                config["network"],
                reference_pno    ,
                wf_fam           ,
                config["fmin"]   ,
            )

            # Loop over the various PN orders to compute the Fisher matrices
            # and errors. 
            for (idx_pno, pno) in enumerate(config["pn_orders"])
                println("\nComputing Fisher matrices for PN order $(pno) with waveform family $(wf_fam)")
                mu = config["mu"][idx_pno]
                sigma = config["sigma"][idx_pno]
                pn_deviation = createSED.createBGRDeviations(catalog, mu, sigma, config["seed"])

                _, _, _, invc, delta_k, dphi_k = createSED.createSEDfromCatalog(
                    catalog          ,
                    config["network"],
                    pno              ,
                    pn_deviation     ,
                    wf_fam           ,
                    config["fmin"]   ,
                    config["seed"]   ;
                    precomputed_snr  = snr,
                    precomputed_isnr = isnr,
                    snr_threshold = config["snr_threshold"]
                )

                write_pn_results_to_hdf5(
                    file             ,
                    pno              ,
                    mu               ,
                    sigma            ,
                    pn_deviation     ,
                    config["network"],
                    wf_fam           ,
                    snr              ,
                    isnr             ,
                    invc             ,
                    delta_k          ,
                    dphi_k           ,
                )
            end
        end
    end

    return output_file
end

#----------------------------------------------------------------------------#
# RUN MAIN 
#----------------------------------------------------------------------------#
function main(args=ARGS)
    if length(args) != 1
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml>"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    config = read_config(config_file)

    return run_fisher_analysis(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
