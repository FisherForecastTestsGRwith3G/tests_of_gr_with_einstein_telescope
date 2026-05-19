using HDF5
using Random
using Dates

include("_config_parser.jl")
include("../scripts_to_run/_hdf5_metadata.jl")
include("../create_single_event_datasets/createSED.jl")

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function _uniform_range(rng::AbstractRNG, bounds::Tuple{Float64, Float64}, n_events::Integer)
    lower, upper = bounds
    return lower .+ (upper - lower) .* rand(rng, n_events)
end

function build_gw150914_like_catalog(config::Dict)
    n_events = config["n_events"]
    event = config["gw150914_like"]
    rng = MersenneTwister(config["seed"])

    spin_range = event["spin_range"]
    theta_cos_range = event["theta_cos_range"]
    phi_range = event["phi_range"]
    iota_cos_range = event["iota_cos_range"]
    psi_range = event["psi_range"]
    t_coal_range = event["t_coal_range"]
    phi_coal_range = event["phi_coal_range"]

    catalog = createSED.BBHCatalog(
        fill(Float64(event["mc"]), n_events),
        fill(Float64(event["eta"]), n_events),
        _uniform_range(rng, spin_range, n_events),
        _uniform_range(rng, spin_range, n_events),
        fill(Float64(event["dL"]), n_events),
        acos.(_uniform_range(rng, theta_cos_range, n_events)),
        _uniform_range(rng, phi_range, n_events),
        acos.(_uniform_range(rng, iota_cos_range, n_events)),
        _uniform_range(rng, psi_range, n_events),
        _uniform_range(rng, t_coal_range, n_events),
        _uniform_range(rng, phi_coal_range, n_events),
        fill(Float64(event["redshift"]), n_events),
    )

    if isfinite(config["z_cut"])
        catalog = createSED.applyRedshiftCut(catalog, config["z_cut"])
        createSED.length(catalog) == n_events ||
            throw(ArgumentError("The configured z_cut removed events from the fixed GW150914-like catalog."))
    end

    return catalog
end

function write_catalog_to_hdf5(file, catalog, config)
    catalog_grp = create_group(file, "bbh_catalog")
    attrs(catalog_grp)["catalog_tag"] = config["catalog_tag"]
    attrs(catalog_grp)["requested_n_events"] = config["n_events"]
    attrs(catalog_grp)["seed"] = config["seed"]
    attrs(catalog_grp)["z_cut"] = config["z_cut"]
    attrs(catalog_grp)["source"] = "GW150914-like fixed-event ensemble"

    param_grp = create_group(catalog_grp, "parameter")
    write(param_grp, "mc", catalog.mc)
    write(param_grp, "eta", catalog.eta)
    write(param_grp, "chi_1", catalog.chi_1)
    write(param_grp, "chi_2", catalog.chi_2)
    write(param_grp, "dL", catalog.dL)
    write(param_grp, "theta", catalog.theta)
    write(param_grp, "phi", catalog.phi)
    write(param_grp, "iota", catalog.iota)
    write(param_grp, "psi", catalog.psi)
    write(param_grp, "t_coal", catalog.t_coal)
    write(param_grp, "phi_coal", catalog.phi_coal)
    write(param_grp, "z", catalog.z)
end

function write_pn_results_to_hdf5(
    file,
    pno,
    mu,
    sigma,
    pn_deviation,
    network,
    wf_fam,
    snr,
    isnr,
    invc,
    delta_k,
    dphi_k,
)
    pnos = createSED.pnoString(pno)
    pno_grp = haskey(file, pnos) ? file[pnos] : create_group(file, pnos)
    attrs(pno_grp)["mu"] = mu
    attrs(pno_grp)["sigma"] = sigma
    if !haskey(pno_grp, "delta_phi_pn")
        write(pno_grp, "delta_phi_pn", pn_deviation)
    end

    network_grp = haskey(pno_grp, network) ? pno_grp[network] : create_group(pno_grp, network)
    wf_grp = haskey(network_grp, wf_fam) ? network_grp[wf_fam] : create_group(network_grp, wf_fam)

    write(wf_grp, "snr", snr)
    write(wf_grp, "isnr", isnr)
    write(wf_grp, "invc", invc)
    write(wf_grp, "delta_k", delta_k)
    write(wf_grp, "dphi_k", dphi_k)
end

function derive_analysis_seed(base_seed::Integer, wf_idx::Integer, pno_idx::Integer, stream_idx::Integer)
    return Int(base_seed + 10_000 * (wf_idx - 1) + 100 * (pno_idx - 1) + stream_idx)
end

function run_fisher_analysis(config::Dict)
    println("Preparing GW150914-like catalog for fmin=$(config["fmin"]) Hz")
    output_file = get_fisher_results_file(config)
    mkpath(dirname(output_file))

    catalog = build_gw150914_like_catalog(config)
    println("Catalog contains $(createSED.length(catalog)) GW150914-like realizations.")
    println("Writing Fisher products to $(output_file)")

    h5open(output_file, "w") do file
        write_top_level_output_metadata!(file)
        attrs(file)["fmin"] = config["fmin"]
        attrs(file)["network"] = config["network"]
        write_catalog_to_hdf5(file, catalog, config)

        reference_pno = first(config["pn_orders"])
        for (idx_wf, wf_fam) in enumerate(config["waveform_families"])
            println("\nProcessing waveform family $(wf_fam) at fmin=$(config["fmin"]) Hz")
            snr, isnr = createSED.computeSNRsFromCatalog(
                catalog,
                config["network"],
                reference_pno,
                wf_fam,
                config["fmin"],
            )

            for (idx_pno, pno) in enumerate(config["pn_orders"])
                println("\nComputing Fisher matrices for PN order $(pno), waveform $(wf_fam)")
                mu = config["mu"][idx_pno]
                sigma = config["sigma"][idx_pno]
                deviation_seed = derive_analysis_seed(config["seed"], 1, idx_pno, 1)
                posterior_seed = derive_analysis_seed(config["seed"], idx_wf, idx_pno, 2)
                pn_deviation = createSED.createBGRDeviations(catalog, mu, sigma, deviation_seed)

                _, _, _, invc, delta_k, dphi_k = createSED.createSEDfromCatalog(
                    catalog,
                    config["network"],
                    pno,
                    pn_deviation,
                    wf_fam,
                    config["fmin"],
                    posterior_seed;
                    precomputed_snr=snr,
                    precomputed_isnr=isnr,
                    snr_threshold=config["snr_threshold"],
                    inspiral_snr_threshold=config["snr_inspiral_threshold"],
                )

                write_pn_results_to_hdf5(
                    file,
                    pno,
                    mu,
                    sigma,
                    pn_deviation,
                    config["network"],
                    wf_fam,
                    snr,
                    isnr,
                    invc,
                    delta_k,
                    dphi_k,
                )
            end
        end
    end

    return output_file
end

function run_all_fmin_analyses(base_config::Dict)
    output_files = String[]
    for config in iter_fmin_configs(base_config)
        push!(output_files, run_fisher_analysis(config))
    end
    return output_files
end

function main(args=ARGS)
    if length(args) != 1
        script_name = basename(@__FILE__)
        throw(ArgumentError("Usage: julia $(script_name) <config_file.toml>"))
    end

    config_file = abspath(args[1])
    println("Using config file: $(config_file)")
    base_config = read_fmin_config(config_file)

    return run_all_fmin_analyses(base_config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
