using TOML
using HDF5
using Printf
using Random
using Dates
using Statistics

include("_config_parser.jl")
include("_hdf5_metadata.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_population_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "population_results_$(config["bootstrap_tag"]).h5")
end

function derive_population_seed(base_seed::Integer, wf_idx::Integer, pno_idx::Integer, stream_idx::Integer)
    return Int(base_seed + 10_000 * (wf_idx - 1) + 100 * (pno_idx - 1) + stream_idx)
end

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#
"""
    write_population_results_hdf5(
        population_results_file,
        config,
        summary_indices,
        pn_indices,
        phi90,
        bootstrap_constraints,
        selected_delta_k,
        selected_dphi_k,
        noninvertible_after_snr,
        bootstrap_observed_fractions,
    )

Write population-analysis results to an HDF5 file.

The file is created at `population_results_file`. Existing files at that path
are overwritten. The function stores results in the following hierarchy:

`/<network>/<waveform_family>/`
- attribute `n_selected`     : number of events selected for the waveform family
- dataset   `summary_indices`: Boolean selection mask for the waveform family

`/<network>/<waveform_family>/<pn_order>/`
- attribute `n_selected`              : number of events selected for this PN order
- attribute `non_invertable_after_snr`: number of events passing SNR cuts but excluded
                                        because the Fisher matrix could not be inverted
- dataset   `pn_indices`              : Boolean selection mask for this PN order
- dataset   `single_event_constraints`: per-event single-event constraints
- dataset   `bootstrap_constraints`   : bootstrap population constraints
- dataset   `bootstrap_observed_fraction`: fraction of the `n_catalog` bootstrap
                                          draw that passed the observation cuts
- dataset   `selected_delta_k`        : selected `delta_k` values used in the population 
                                        analysis.
- dataset   `selected_dphi_k`         : selected `dphi_k` values used in the population
                                        analysis

Here `<network>` is `config["network"]`, each `<waveform_family>` comes from
`config["waveform_families"]`, and each `<pn_order>` is the string produced by
`createSED.pnoString(pno)` for `pno` in `config["pn_orders"]`.

Arguments
---------
- `population_results_file`: Output HDF5 file path.
- `config`                 : Analysis configuration dictionary containing the network, 
                             waveform families, PN orders, and output settings.
- `summary_indices`        : Per-waveform-family event masks after all selection cuts.
- `pn_indices`             : Per-waveform-family, per-PN-order event masks.
- `phi90`                  : Per-waveform-family, per-PN-order vectors of single-event 
                             constraints.
- `bootstrap_constraints`  : Per-waveform-family, per-PN-order vectors of bootstrap 
                             constraints.
- `selected_delta_k`       : Per-waveform-family, per-PN-order vectors of selected `delta_k` 
                             values.
- `selected_dphi_k`        : Per-waveform-family, per-PN-order vectors of selected `dphi_k` 
                             values.
- `noninvertible_after_snr`: Per-waveform-family, per-PN-order masks for events that 
                            pass SNR cuts but fail the invertibility requirement.
- `bootstrap_observed_fractions`: Per-waveform-family, per-PN-order vectors containing
                                  the observed-event fraction for each bootstrap sample.
"""
function write_population_results_hdf5(
    population_results_file::AbstractString,
    config::Dict,
    summary_indices::Dict{String, BitVector},
    pn_indices::Dict{String, Dict{String, BitVector}},
    phi90::Dict{String, Dict{String, Vector{Float64}}},
    bootstrap_constraints::Dict{String, Dict{String, Vector{Float64}}},
    selected_delta_k::Dict{String, Dict{String, Vector{Float64}}},
    selected_dphi_k::Dict{String, Dict{String, Vector{Float64}}},
    noninvertible_after_snr::Dict{String, Dict{String, BitVector}},
    bootstrap_observed_fractions::Dict{String, Dict{String, Vector{Float64}}},
)
    mkpath(dirname(population_results_file))

    h5open(population_results_file, "w") do file
        write_top_level_output_metadata!(file)
        attrs(file)["select_before_bootstrap"] = config["select_before_bootstrap"]
        attrs(file)["bootstrap_selection_mode"] = config["select_before_bootstrap"] ? "select_before_bootstrap" : "bootstrap_before_select"
        network_group = create_group(file, config["network"])

        for wf_fam in config["waveform_families"]
            wf_grp = create_group(network_group, wf_fam)
            attrs(wf_grp)["n_selected"] = sum(summary_indices[wf_fam])
            write(wf_grp, "summary_indices", Vector{Bool}(summary_indices[wf_fam]))

            for pno in config["pn_orders"]
                pno_grp = create_group(wf_grp, createSED.pnoString(pno))
                attrs(pno_grp)["n_selected"] = sum(pn_indices[wf_fam][pno])
                attrs(pno_grp)["non_invertable_after_snr"] = sum(noninvertible_after_snr[wf_fam][pno])

                write(pno_grp, "pn_indices", Vector{Bool}(pn_indices[wf_fam][pno]))
                write(pno_grp, "single_event_constraints", phi90[wf_fam][pno])
                write(pno_grp, "bootstrap_constraints", bootstrap_constraints[wf_fam][pno])
                write(pno_grp, "bootstrap_observed_fraction", bootstrap_observed_fractions[wf_fam][pno])
                write(pno_grp, "selected_delta_k", selected_delta_k[wf_fam][pno])
                write(pno_grp, "selected_dphi_k", selected_dphi_k[wf_fam][pno])
            end
        end
    end
end

"""
    print_population_workflow_header(mode_label, config)

Print the analysis mode and the common selection/bootstrap settings.
"""
function print_population_workflow_header(mode_label::AbstractString, config::Dict)
    println("\nPopulation analysis workflow: $(mode_label)")
    println("Applied SNR cuts: snr > $(config["snr_threshold"]), isnr > $(config["snr_inspiral_threshold"])")
    println("Bootstrap settings: n_catalog = $(config["n_catalog"]), n_sample = $(config["n_sample"])")
end

"""
    print_select_before_bootstrap_summary(...)

Print a per-waveform, per-PN-order summary for the workflow that applies the
selection before the bootstrap stage.
"""
function print_select_before_bootstrap_summary(
    wf_fam::AbstractString,
    pno::AbstractString,
    total_events::Integer,
    pn_selected::Integer,
    shared_selected::Integer,
    noninvertible_count::Integer,
    n_constraints::Integer,
)
    println("\nWaveform family: $(wf_fam), PN order: $(pno)")
    println("  Total events loaded: $(total_events)")
    println("  Events passing SNR and invertibility cuts for this PN order: $(pn_selected)")
    println("  SNR-selected but non-invertible events: $(noninvertible_count)")
    println("  Events retained after cross-PN selection: $(shared_selected)")
    println(" *Running phi90 and bootstrap on the shared selected sample.")
    println("  Bootstrap constraints generated: $(n_constraints)")
end

"""
    print_bootstrap_before_select_summary(...)

Print a per-waveform, per-PN-order summary for the workflow that bootstraps the
catalog before applying the observation selection.
"""
function print_bootstrap_before_select_summary(
    wf_fam::AbstractString,
    pno::AbstractString,
    total_events::Integer,
    pn_selected::Integer,
    shared_selected::Integer,
    noninvertible_count::Integer,
    observed_fractions::Vector{Float64},
)
    empty_draws = count(iszero, observed_fractions)
    observed_fraction_ci90 = quantile(observed_fractions, [0.05, 0.95])
    observed_fraction_median = median(observed_fractions)
    observed_fraction_mean = mean(observed_fractions)
    observed_fraction_std = std(observed_fractions)
    format_percent(value) = @sprintf("%.2f%%", 100 * value)
    println("\nWaveform family: $(wf_fam), PN order: $(pno)")
    println("  Total events loaded: $(total_events)")
    println("  Events passing SNR and invertibility cuts for this PN order: $(pn_selected)")
    println("  SNR-selected but non-invertible events: $(noninvertible_count)")
    println("  Events retained after cross-PN selection for phi90: $(shared_selected)")
    println(" *Bootstrap draws were taken from the full catalog and filtered afterwards.")
    println("  Observed fraction per draw (median -/+ lower/upper 90% CI errr): $(format_percent(observed_fraction_median)) - $(format_percent(observed_fraction_median - observed_fraction_ci90[1])) + $(format_percent(observed_fraction_ci90[2] - observed_fraction_median))")
    println("  Observed fraction per draw (mean +- std): $(format_percent(observed_fraction_mean)) +- $(format_percent(observed_fraction_std))")
    println("  Observed fraction per draw ([0.05,0.5,0.95]-quantile): [$(format_percent(observed_fraction_ci90[1])), $(format_percent(observed_fraction_median)), $(format_percent(observed_fraction_ci90[2]))]")
    println("  Empty observed draws: $(empty_draws) / $(length(observed_fractions))")
end

"""
Calculates the 90%Credible upper bound on the GR Deviation in the limit where 
the deviation distribution is not the hierarchical one, but simple the one obtained
by multiplying all the likelihoods. Physically, this amounts to the assumption
that all GR deviation manifest in exactly the same way. 

The constrained obtained in this way will always be tighter than the hierarchical 
constrained, with the gaussian population model. 
"""
function prodLPopDeviationConstraint(dphi_k::Vector{Float64}, delta_k::Vector{Float64})

    mu_eff, sig_eff = HierDist.muStdEff4ProdNormal(dphi_k, delta_k)
    return HierDist.estimate0Symmetric90CiGaussian(mu_eff, sig_eff)
end 

"""
Estimate a population-level GR-deviation constraint from bootstrap resamples of
single-event Gaussian likelihood summaries.

For each bootstrap realization, the function draws `n_catalog` distinct event
indices without replacement from the input arrays, evaluates the requested
population-constraint method on that subset, and stores the resulting
constraint value. Repeating this `n_samples` times returns a bootstrap sample
of population constraints.

Population Constraint Methods
----------------------------
- "prodL"    : The deviation distribution is calculated as the product of 
               individual likelihoods of single event constraints
- "Hier"     : Hierarchical combination of constraints, assuming a gaussian
               Population model for the deviations. 

Arguments
---------
- `dphi_k`   : Mean values of the single-event Gaussian constraints.
- `delta_k`  : Standard deviations of the single-event Gaussian constraints.
- `n_catalog`: Number of events drawn per bootstrap realization.
- `n_samples`: Number of bootstrap realizations.
- `method`   : Population-constraint method. Currently supported: `"prodL"`.
- `seed`     : Random seed used for the bootstrap draws.

Returns
-------
- `Vector{Float64}`: One bootstrap population constraint per realization.
"""
function bootstrapPopConstraint(
    dphi_k::Vector{Float64}, 
    delta_k::Vector{Float64}, 
    n_catalog::Int, 
    n_samples::Int, 
    method::String;
    seed::Int=1234,
    )

    length(dphi_k) == length(delta_k) ||
        throw(ArgumentError("`dphi_k` and `delta_k` must have the same length."))
    n_catalog > 0 || throw(ArgumentError("`n_catalog` must be positive."))
    n_samples > 0 || throw(ArgumentError("`n_samples` must be positive."))
    n_catalog <= length(dphi_k) ||
        throw(ArgumentError("`n_catalog` cannot exceed the number of available events."))

    n_combinations = binomial(big(length(dphi_k)), big(n_catalog))
    n_combinations >= 1000 ||
        throw(ArgumentError(
            "Too few distinct event combinations for bootstrap: choose($(length(dphi_k)), $(n_catalog)) = $(n_combinations) < 1000."
        ))

    if method == "prodL"
        constraint(mu, sig) = prodLPopDeviationConstraint(mu, sig)
    else
        throw(ArgumentError("Unknown bootstrap method: $(method)"))
    end

    rng = MersenneTwister(seed)
    constraints = Vector{Float64}(undef, n_samples)
    n_events = length(dphi_k)

    for idx_sample in 1:n_samples
        selected_indices = randperm(rng, n_events)[1:n_catalog]
        constraints[idx_sample] = constraint(
            dphi_k[selected_indices],
            delta_k[selected_indices],
        )
    end

    return constraints
end

"""
Estimate population-level GR-deviation constraints when bootstrap resampling is
performed before applying the observation selection.

For each bootstrap realization, the function draws `n_catalog` distinct catalog
events without replacement, restricts that draw to the entries marked `true` in
`selection_index`, and evaluates the requested population-constraint method on
the surviving observed subset. It returns both the resulting bootstrap
constraints and the observed fraction in each draw, defined as the number of
selected events divided by `n_catalog`.

If a bootstrap draw contains no observed events, the corresponding constraint
is stored as `NaN`.

Arguments
---------
- `dphi_k`         : Mean values of the single-event Gaussian constraints for the full catalog.
- `delta_k`        : Standard deviations of the single-event Gaussian constraints.
- `selection_index`: Boolean mask indicating which catalog events satisfy the observation cuts.
- `n_catalog`      : Number of catalog events drawn per bootstrap realization.
- `n_samples`      : Number of bootstrap realizations.
- `method`         : Population-constraint method. Currently supported: `"prodL"`.
- `seed`           : Random seed used for the bootstrap draws.

Returns
-------
- `(constraints, observed_fractions)`: Tuple of vectors containing the bootstrap
  constraint values and the observed-event fraction for each bootstrap sample.
"""
function bootstrapPopConstraintAfterSelection(
    dphi_k::Vector{Float64},
    delta_k::Vector{Float64},
    selection_index::BitVector,
    n_catalog::Int,
    n_samples::Int,
    method::String;
    seed::Int=1234,
)
    length(dphi_k) == length(delta_k) ||
        throw(ArgumentError("`dphi_k` and `delta_k` must have the same length."))
    length(selection_index) == length(dphi_k) ||
        throw(ArgumentError("`selection_index` must have the same length as `dphi_k` and `delta_k`."))
    n_catalog > 0 || throw(ArgumentError("`n_catalog` must be positive."))
    n_samples > 0 || throw(ArgumentError("`n_samples` must be positive."))
    n_catalog <= length(dphi_k) ||
        throw(ArgumentError("`n_catalog` cannot exceed the number of available catalog events."))

    if method == "prodL"
        constraint(mu, sig) = prodLPopDeviationConstraint(mu, sig)
    else
        throw(ArgumentError("Unknown bootstrap method: $(method)"))
    end

    rng = MersenneTwister(seed)
    constraints = Vector{Float64}(undef, n_samples)
    observed_fractions = Vector{Float64}(undef, n_samples)
    n_events = length(dphi_k)
    for idx_sample in 1:n_samples
        sampled_indices = randperm(rng, n_events)[1:n_catalog]
        observed_indices = sampled_indices[selection_index[sampled_indices]]
        observed_fractions[idx_sample] = length(observed_indices) / n_catalog

        if isempty(observed_indices)
            constraints[idx_sample] = NaN
        else
            constraints[idx_sample] = constraint(
                dphi_k[observed_indices],
                delta_k[observed_indices],
            )
        end
    end

    return constraints, observed_fractions
end


"""
Read the Fisher-analysis products needed for the population stage and build the
derived selection masks and event-level summaries.

For each waveform family and PN order, the function loads `snr`, `isnr`,
`invc`, `delta_k`, and `dphi_k` from the Fisher-results file. It constructs:

- `pn_indices`: events that pass the SNR, inspiral-SNR, and invertibility cuts
  for each PN order
- `summary_indices`: events that satisfy those cuts across all configured PN
  orders for a waveform family
- `noninvertible_after_snr`: events that pass the SNR cuts but fail the
  invertibility requirement

It also loads and returns the full per-event `delta_k` and `dphi_k` arrays for
each waveform family and PN order so that downstream workflows can either
apply the selection before bootstrapping or after bootstrapping.

Arguments
---------
- `file`  : Open HDF5 handle to the Fisher-results file.
- `config`: Analysis configuration containing the network, waveform families,
            PN orders, and selection thresholds.

Returns
-------
- `(pn_indices, summary_indices, noninvertible_after_snr, delta_k_data, dphi_k_data)`:
  tuple of dictionaries containing the selection masks and full event-level
  summaries used by the population-analysis workflows.
"""
function collect_population_inputs(file, config::Dict)
    pn_indices = Dict{String, Dict{String, BitVector}}()
    summary_indices = Dict{String, BitVector}()
    noninvertible_after_snr = Dict{String, Dict{String, BitVector}}()
    delta_k_data = Dict{String, Dict{String, Vector{Float64}}}()
    dphi_k_data = Dict{String, Dict{String, Vector{Float64}}}()

    for wf_fam in config["waveform_families"]
        summary_index = nothing

        for pno in config["pn_orders"]
            pno_grp = createSED.pnoString(pno)
            wf_grp = file[pno_grp][config["network"]][wf_fam]

            snr  = read(wf_grp, "snr")
            isnr = read(wf_grp, "isnr")
            invc = read(wf_grp, "invc")

            selection_index = BitVector(
                (snr .> config["snr_threshold"]) .&
                (isnr .> config["snr_inspiral_threshold"]) .&
                invc
            )
            snr_not_invc_index = BitVector(
                (snr .> config["snr_threshold"]) .&
                (isnr .> config["snr_inspiral_threshold"]) .&
                .!invc
            )

            if !haskey(pn_indices, wf_fam)
                pn_indices[wf_fam] = Dict{String, BitVector}()
            end
            pn_indices[wf_fam][pno] = selection_index

            if !haskey(noninvertible_after_snr, wf_fam)
                noninvertible_after_snr[wf_fam] = Dict{String, BitVector}()
            end
            noninvertible_after_snr[wf_fam][pno] = snr_not_invc_index

            if isnothing(summary_index)
                summary_index = copy(selection_index)
            else
                summary_index .&= selection_index
            end
        end

        summary_indices[wf_fam] = summary_index
    end

    for wf_fam in config["waveform_families"]
        delta_k_data[wf_fam] = Dict{String, Vector{Float64}}()
        dphi_k_data[wf_fam] = Dict{String, Vector{Float64}}()
        for pno in config["pn_orders"]
            pno_grp = createSED.pnoString(pno)
            wf_grp = file[pno_grp][config["network"]][wf_fam]

            delta_k_data[wf_fam][pno] = read(wf_grp, "delta_k")
            dphi_k_data[wf_fam][pno] = read(wf_grp, "dphi_k")
        end
    end

    return pn_indices, summary_indices, noninvertible_after_snr, delta_k_data, dphi_k_data
end

"""
Run the population-analysis workflow that applies the event-selection mask
before bootstrapping the catalog-level constraint.

The function first loads the shared population-analysis inputs with
`collect_population_inputs`, then restricts `delta_k` and `dphi_k` to the
events that pass the combined selection across all configured PN orders for
each waveform family. It computes the direct 90% interval estimate for the
selected events and evaluates bootstrap population constraints on those
already-selected samples.

Arguments
---------
- `file`  : Open HDF5 handle to the Fisher-results file.
- `config`: Analysis configuration containing waveform families, PN orders,
            bootstrap settings, and selection thresholds.

Returns
-------
- Named tuple containing the shared selection masks, the selected event-level
  `delta_k` and `dphi_k` samples, per-configuration `phi90` estimates,
  bootstrap constraints, and the bootstrap observed fractions.
"""
function run_population_analysis_select_before_bootstrap(file, config::Dict; verbose::Bool=true)
    pn_indices, summary_indices, noninvertible_after_snr, delta_k_data, dphi_k_data =
        collect_population_inputs(file, config)

    selected_delta_k = Dict{String, Dict{String, Vector{Float64}}}()
    selected_dphi_k = Dict{String, Dict{String, Vector{Float64}}}()
    phi90 = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_constraints = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_observed_fractions = Dict{String, Dict{String, Vector{Float64}}}()

    verbose && print_population_workflow_header("select_before_bootstrap", config)

    for (idx_wf, wf_fam) in enumerate(config["waveform_families"])
        selected_delta_k[wf_fam] = Dict{String, Vector{Float64}}()
        selected_dphi_k[wf_fam] = Dict{String, Vector{Float64}}()
        phi90[wf_fam] = Dict{String, Vector{Float64}}()
        bootstrap_constraints[wf_fam] = Dict{String, Vector{Float64}}()
        bootstrap_observed_fractions[wf_fam] = Dict{String, Vector{Float64}}()
        selection_index = summary_indices[wf_fam]

        for (idx_pno, pno) in enumerate(config["pn_orders"])
            selected_delta_k[wf_fam][pno] = delta_k_data[wf_fam][pno][selection_index]
            selected_dphi_k[wf_fam][pno] = dphi_k_data[wf_fam][pno][selection_index]

            phi90[wf_fam][pno] = HierDist.estimate0Symmetric90CiGaussian(
                selected_dphi_k[wf_fam][pno],
                selected_delta_k[wf_fam][pno]
            )

            bootstrap_constraints[wf_fam][pno] = bootstrapPopConstraint(
                selected_dphi_k[wf_fam][pno],
                selected_delta_k[wf_fam][pno],
                config["n_catalog"],
                config["n_sample"],
                "prodL",
                seed=derive_population_seed(config["seed"], idx_wf, idx_pno, 1),
            )
            bootstrap_observed_fractions[wf_fam][pno] = ones(Float64, config["n_sample"])

            if verbose
                print_select_before_bootstrap_summary(
                    wf_fam,
                    pno,
                    length(delta_k_data[wf_fam][pno]),
                    sum(pn_indices[wf_fam][pno]),
                    length(selected_delta_k[wf_fam][pno]),
                    sum(noninvertible_after_snr[wf_fam][pno]),
                    length(bootstrap_constraints[wf_fam][pno]),
                )
            end
        end
    end

    return (
        summary_indices=summary_indices,
        pn_indices=pn_indices,
        phi90=phi90,
        bootstrap_constraints=bootstrap_constraints,
        selected_delta_k=selected_delta_k,
        selected_dphi_k=selected_dphi_k,
        noninvertible_after_snr=noninvertible_after_snr,
        bootstrap_observed_fractions=bootstrap_observed_fractions,
    )
end

"""
Run the population-analysis workflow that bootstraps the catalog-level
constraint before applying the final cross-PN event selection.

The function first loads the shared population-analysis inputs with
`collect_population_inputs`. For each waveform family and PN order, it keeps
the PN-specific selection mask for the bootstrap stage, computes bootstrap
constraints on the full event set with selection applied inside the resampling
step, and separately stores the `delta_k` and `dphi_k` samples that satisfy
the combined selection across all configured PN orders for direct `phi90`
estimation.

Arguments
---------
- `file`  : Open HDF5 handle to the Fisher-results file.
- `config`: Analysis configuration containing waveform families, PN orders,
            bootstrap settings, random seed, and selection thresholds.

Returns
-------
- Named tuple containing the shared selection masks, the selected event-level
  `delta_k` and `dphi_k` samples, per-configuration `phi90` estimates,
  bootstrap constraints, and the bootstrap observed fractions.
"""
function run_population_analysis_bootstrap_before_select(file, config::Dict; verbose::Bool=true)
    pn_indices, summary_indices, noninvertible_after_snr, delta_k_data, dphi_k_data =
        collect_population_inputs(file, config)

    selected_delta_k = Dict{String, Dict{String, Vector{Float64}}}()
    selected_dphi_k = Dict{String, Dict{String, Vector{Float64}}}()
    phi90 = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_constraints = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_observed_fractions = Dict{String, Dict{String, Vector{Float64}}}()

    verbose && print_population_workflow_header("bootstrap_before_select", config)

    for (idx_wf, wf_fam) in enumerate(config["waveform_families"])
        selected_delta_k[wf_fam] = Dict{String, Vector{Float64}}()
        selected_dphi_k[wf_fam] = Dict{String, Vector{Float64}}()
        phi90[wf_fam] = Dict{String, Vector{Float64}}()
        bootstrap_constraints[wf_fam] = Dict{String, Vector{Float64}}()
        bootstrap_observed_fractions[wf_fam] = Dict{String, Vector{Float64}}()
        summary_index = summary_indices[wf_fam]

        for (idx_pno, pno) in enumerate(config["pn_orders"])
            delta_k = delta_k_data[wf_fam][pno]
            dphi_k = dphi_k_data[wf_fam][pno]
            selection_index = pn_indices[wf_fam][pno]

            selected_delta_k[wf_fam][pno] = delta_k[summary_index]
            selected_dphi_k[wf_fam][pno] = dphi_k[summary_index]

            phi90[wf_fam][pno] = HierDist.estimate0Symmetric90CiGaussian(
                selected_dphi_k[wf_fam][pno],
                selected_delta_k[wf_fam][pno]
            )

            constraints, observed_fractions = bootstrapPopConstraintAfterSelection(
                dphi_k,
                delta_k,
                selection_index,
                config["n_catalog"],
                config["n_sample"],
                "prodL";
                seed=derive_population_seed(config["seed"], idx_wf, idx_pno, 2),
            )
            bootstrap_constraints[wf_fam][pno] = constraints
            bootstrap_observed_fractions[wf_fam][pno] = observed_fractions

            if verbose
                print_bootstrap_before_select_summary(
                    wf_fam,
                    pno,
                    length(delta_k),
                    sum(selection_index),
                    length(selected_delta_k[wf_fam][pno]),
                    sum(noninvertible_after_snr[wf_fam][pno]),
                    observed_fractions,
                )
            end
        end
    end

    return (
        summary_indices=summary_indices,
        pn_indices=pn_indices,
        phi90=phi90,
        bootstrap_constraints=bootstrap_constraints,
        selected_delta_k=selected_delta_k,
        selected_dphi_k=selected_dphi_k,
        noninvertible_after_snr=noninvertible_after_snr,
        bootstrap_observed_fractions=bootstrap_observed_fractions,
    )
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
function run_population_analysis(config::Dict; verbose::Bool=true)
    fisher_results_file = get_fisher_results_file(config)
    population_results_file = get_population_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    println("Using Fisher results file: $(fisher_results_file)")

    h5open(fisher_results_file, "r") do file
        results = if config["select_before_bootstrap"]
            run_population_analysis_select_before_bootstrap(file, config) 
        else
            run_population_analysis_bootstrap_before_select(file, config)
        end
        
        write_population_results_hdf5(
            population_results_file,    
            config,
            results.summary_indices,
            results.pn_indices,
            results.phi90,
            results.bootstrap_constraints,
            results.selected_delta_k,
            results.selected_dphi_k,
            results.noninvertible_after_snr,
            results.bootstrap_observed_fractions,
        )
    end

    println("Population results written to $(population_results_file)")
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

    return run_population_analysis(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
