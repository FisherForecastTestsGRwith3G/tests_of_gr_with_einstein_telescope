using TOML

function _require_section(raw_config::Dict, section::AbstractString)
    haskey(raw_config, section) ||
        throw(ArgumentError("Missing required config section: [$(section)]"))
    return raw_config[section]
end

function _parse_range(section::Dict, key::AbstractString, default::Tuple{<:Real, <:Real})
    values = haskey(section, key) ? Float64.(section[key]) : collect(default)
    length(values) == 2 ||
        throw(ArgumentError("`$(key)` must contain exactly two values."))
    values[1] <= values[2] ||
        throw(ArgumentError("`$(key)` lower bound must not exceed the upper bound."))
    return (values[1], values[2])
end

function _parse_optional_limits(section::Dict, key::AbstractString)
    haskey(section, key) || return nothing
    values = Float64.(section[key])
    length(values) == 2 ||
        throw(ArgumentError("`plots.$(key)` must contain exactly two values."))
    all(values .> 0.0) ||
        throw(ArgumentError("`plots.$(key)` must be positive because the plot uses log axes."))
    return (values[1], values[2])
end

function format_fmin_label(fmin::Real)
    value = Float64(fmin)
    if isapprox(value, round(value); atol=1e-9)
        return string(round(Int, value), "Hz")
    end
    return replace(string(value), "." => "p") * "Hz"
end

function read_fmin_config(config_file::AbstractString)
    config_file = abspath(config_file)
    isfile(config_file) || throw(ArgumentError("Config file does not exist: $(config_file)"))

    raw_config = TOML.parsefile(config_file)

    general = _require_section(raw_config, "general")
    catalog = _require_section(raw_config, "catalog")
    event = _require_section(raw_config, "gw150914_like")
    deviations = _require_section(raw_config, "deviations")
    detectors = _require_section(raw_config, "detectors")
    fisher = _require_section(raw_config, "fisher")
    population = _require_section(raw_config, "population")
    fmin_comparison = _require_section(raw_config, "fmin_comparison")
    bootstrap = get(raw_config, "bootstrap", Dict{String, Any}())
    plots = get(raw_config, "plots", Dict{String, Any}())

    config = Dict{String, Any}()
    config["config_file"] = config_file
    config["seed"] = Int(general["seed"])
    config["outdir"] = String(get(general, "outdir", get(catalog, "outdir", "./results/data/")))

    config["catalog_tag"] = String(catalog["catalog_tag"])
    config["n_events"] = Int(catalog["n_events"])
    config["z_cut"] = haskey(catalog, "z_cut") ? Float64(catalog["z_cut"]) : Inf

    config["gw150914_like"] = Dict{String, Any}(
        "mc" => Float64(get(event, "mc", 30.0)),
        "eta" => Float64(get(event, "eta", 0.247)),
        "dL" => Float64(get(event, "dL", 0.4)),
        "redshift" => Float64(get(event, "redshift", 0.09)),
        "spin_range" => _parse_range(event, "spin_range", (-0.3, 0.3)),
        "theta_cos_range" => _parse_range(event, "theta_cos_range", (-1.0, 1.0)),
        "phi_range" => _parse_range(event, "phi_range", (0.0, 2 * pi)),
        "iota_cos_range" => _parse_range(event, "iota_cos_range", (-1.0, cosd(135.0))),
        "psi_range" => _parse_range(event, "psi_range", (0.0, pi)),
        "t_coal_range" => _parse_range(event, "t_coal_range", (0.0, 1.0)),
        "phi_coal_range" => _parse_range(event, "phi_coal_range", (0.0, 2 * pi)),
    )

    config["pn_orders"] = String.(deviations["pn_orders"])
    config["mu"] = Float64.(deviations["mu"])
    config["sigma"] = Float64.(deviations["sigma"])
    length(config["pn_orders"]) == length(config["mu"]) ||
        throw(ArgumentError("`deviations.mu` must have the same length as `deviations.pn_orders`."))
    length(config["pn_orders"]) == length(config["sigma"]) ||
        throw(ArgumentError("`deviations.sigma` must have the same length as `deviations.pn_orders`."))

    network_config = detectors["network"]
    network_config isa AbstractString ||
        throw(ArgumentError("`detectors.network` must be a single string, e.g. network = \"ETS\"."))
    config["network"] = String(network_config)

    waveform_config = fisher["waveform"]
    waveform_config isa AbstractVector ||
        throw(ArgumentError("`fisher.waveform` must be a list, e.g. waveform = [\"PhenomHM\"]."))
    config["waveform_families"] = String.(waveform_config)

    config["fmin_values"] = Float64.(fmin_comparison["fmin_values"])
    isempty(config["fmin_values"]) &&
        throw(ArgumentError("`fmin_comparison.fmin_values` must contain at least one value."))
    config["fmin"] = first(config["fmin_values"])

    config["snr_threshold"] = Float64(population["snr_threshold"])
    config["snr_inspiral_threshold"] = Float64(get(population, "snr_inspiral_threshold", 0.0))
    config["select_before_bootstrap"] = Bool(get(population, "select_before_bootstrap", true))

    config["bootstrap_outdir"] = String(get(bootstrap, "outdir", config["outdir"]))
    config["bootstrap_tag"] = String(get(bootstrap, "bootstrap_tag", config["catalog_tag"]))
    config["n_catalog"] = Int(get(bootstrap, "n_catalog", 1))
    config["n_sample"] = Int(get(bootstrap, "n_sample", 1000))

    config["plot_outdir"] = String(get(plots, "outdir", "./results/plots/"))
    config["plot_tag"] = String(get(plots, "plot_tag", config["catalog_tag"]))
    config["plot_waveform"] = String(get(plots, "waveform", first(config["waveform_families"])))
    config["plot_events_per_realization"] = Int(get(plots, "events_per_realization", 1))
    config["plot_compute_log_summary"] = Bool(get(plots, "compute_mean_std_dev_in_log_space", true))
    config["plot_y_limits_minus_one"] = _parse_optional_limits(plots, "y_limits_minus_one")
    config["plot_y_limits_remaining"] = _parse_optional_limits(plots, "y_limits_remaining")

    config["plot_waveform"] in config["waveform_families"] ||
        throw(ArgumentError("`plots.waveform` must be one of `fisher.waveform`."))
    config["plot_events_per_realization"] > 0 ||
        throw(ArgumentError("`plots.events_per_realization` must be positive."))

    return config
end

function config_for_fmin(base_config::Dict, fmin::Real)
    config = deepcopy(base_config)
    label = format_fmin_label(fmin)
    config["fmin"] = Float64(fmin)
    config["fmin_label"] = label
    config["catalog_tag"] = "$(base_config["catalog_tag"])_fmin_$(label)"
    config["bootstrap_tag"] = "$(base_config["bootstrap_tag"])_fmin_$(label)"
    return config
end

function iter_fmin_configs(base_config::Dict)
    return [config_for_fmin(base_config, fmin) for fmin in base_config["fmin_values"]]
end
