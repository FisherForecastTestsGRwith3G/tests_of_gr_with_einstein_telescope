function read_config(config_file::AbstractString)
    config_file = abspath(config_file)
    isfile(config_file) || throw(ArgumentError("Config file does not exist: $(config_file)"))

    raw_config = TOML.parsefile(config_file)

    required_sections = ("general", "catalog", "deviations", "detectors", "fisher", "population", "bootstrap")
    for section in required_sections
        haskey(raw_config, section) || throw(ArgumentError("Missing required config section: [$(section)]"))
    end

    config = Dict{String, Any}()

    config["config_file"] = config_file
    config["seed"] = Int(raw_config["general"]["seed"])
    if haskey(raw_config["general"], "outdir")
        config["outdir"] = String(raw_config["general"]["outdir"])
    elseif haskey(raw_config["catalog"], "outdir")
        config["outdir"] = String(raw_config["catalog"]["outdir"])
    else
        config["outdir"] = "."
    end

    config["catalog_tag"] = String(raw_config["catalog"]["catalog_tag"])
    config["catalog_outdir"] = haskey(raw_config["catalog"], "outdir") ? String(raw_config["catalog"]["outdir"]) : config["outdir"]
    config["n_events"] = Int(raw_config["catalog"]["n_events"])
    config["z_cut"] = haskey(raw_config["catalog"], "z_cut") ? Float64(raw_config["catalog"]["z_cut"]) : Inf

    config["pn_orders"] = String.(raw_config["deviations"]["pn_orders"])
    config["mu"] = Float64.(raw_config["deviations"]["mu"])
    config["sigma"] = Float64.(raw_config["deviations"]["sigma"])

    length(config["pn_orders"]) == length(config["mu"]) ||
        throw(ArgumentError("`deviations.mu` must have the same length as `deviations.pn_orders`."))
    length(config["pn_orders"]) == length(config["sigma"]) ||
        throw(ArgumentError("`deviations.sigma` must have the same length as `deviations.pn_orders`."))

    network_config = raw_config["detectors"]["network"]
    network_config isa AbstractString ||
        throw(ArgumentError("`detectors.network` must be given as a single string, e.g. network = \"LHV\""))
    config["network"] = String(network_config)

    waveform_config = raw_config["fisher"]["waveform"]
    waveform_config isa AbstractVector ||
        throw(ArgumentError("`fisher.waveform` must be given as a list, e.g. waveform = [\"PhenomHM\"]"))
    config["waveform_families"] = String.(waveform_config)
    config["fmin"] = Float64(raw_config["fisher"]["fmin"])
    config["snr_threshold"] = Float64(raw_config["population"]["snr_threshold"])
    config["snr_inspiral_threshold"] = Float64(raw_config["population"]["snr_inspiral_threshold"])
    config["select_before_bootstrap"] = haskey(raw_config["population"], "select_before_bootstrap") ?
        Bool(raw_config["population"]["select_before_bootstrap"]) : true

    config["bootstrap_outdir"] = haskey(raw_config["bootstrap"], "outdir") ? String(raw_config["bootstrap"]["outdir"]) : config["outdir"]
    config["bootstrap_tag"] = haskey(raw_config["bootstrap"], "bootstrap_tag") ? String(raw_config["bootstrap"]["bootstrap_tag"]) : config["catalog_tag"]
    config["n_catalog"] = Int(raw_config["bootstrap"]["n_catalog"])
    config["n_sample"] = Int(raw_config["bootstrap"]["n_sample"])

    if haskey(raw_config, "scaling")
        scaling = raw_config["scaling"]
        config["scaling_outdir"] = haskey(scaling, "outdir") ? String(scaling["outdir"]) : config["outdir"]
        config["scaling_tag"] = haskey(scaling, "scaling_tag") ? String(scaling["scaling_tag"]) : config["bootstrap_tag"]
        config["scaling_n_max"] = haskey(scaling, "n_max") ? Int(scaling["n_max"]) : config["n_catalog"]
        config["scaling_n_steps"] = haskey(scaling, "n_steps") ? Int(scaling["n_steps"]) : 10
        config["scaling_n_sample"] = haskey(scaling, "n_sample") ? Int(scaling["n_sample"]) : config["n_sample"]
    else
        config["scaling_outdir"] = config["outdir"]
        config["scaling_tag"] = config["bootstrap_tag"]
        config["scaling_n_max"] = config["n_catalog"]
        config["scaling_n_steps"] = 10
        config["scaling_n_sample"] = config["n_sample"]
        config["scaling_method"] = "prodL"
    end

    if haskey(raw_config, "plots")
        config["plot_outdir"] = haskey(raw_config["plots"], "outdir") ? String(raw_config["plots"]["outdir"]) : config["outdir"]
        config["plot_tag"] = haskey(raw_config["plots"], "plot_tag") ? String(raw_config["plots"]["plot_tag"]) : config["catalog_tag"]
    else
        config["plot_outdir"] = config["outdir"]
        config["plot_tag"] = config["catalog_tag"]
    end

    return config
end
