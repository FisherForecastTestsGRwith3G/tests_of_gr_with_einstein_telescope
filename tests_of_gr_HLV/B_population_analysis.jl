using TOML
using HDF5
using Random

include("_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_population_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "population_results_$(config["bootstrap_tag"]).h5")
end

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
)
    mkpath(dirname(population_results_file))

    h5open(population_results_file, "w") do file
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
                write(pno_grp, "selected_delta_k", selected_delta_k[wf_fam][pno])
                write(pno_grp, "selected_dphi_k", selected_dphi_k[wf_fam][pno])
            end
        end
    end
end

"""
    print_selection_summary(config, pn_indices, summary_indices, noninvertible_after_snr)

Print a summary of event selection counts for each waveform family and PN order.

The summary includes the number of selected events after PN-order-specific cuts,
the number of events that pass the SNR cuts but become non-invertible, and the
final number of events retained per waveform family. The function also computes
and returns the event-selection mask shared across all waveform families.

Arguments
---------
- `config`                 : Analysis configuration containing waveform families, 
                             PN orders, and SNR thresholds.
- `pn_indices`             : Per-waveform-family and per-PN-order selection masks.
- `summary_indices`        : Final per-waveform-family selection masks.
- `noninvertible_after_snr`: Per-waveform-family and per-PN-order masks for
                             events that pass SNR cuts but are non-invertible.

Returns
-------
- `shared_index`           : Bit vector containing the events selected in every 
                             waveform family.
"""
function print_selection_summary(config::Dict, pn_indices::Dict{String, Dict{String, BitVector}}, summary_indices::Dict{String, BitVector}, noninvertible_after_snr::Dict{String, Dict{String, BitVector}})
    println("\nSelection summary")
    println("Applied SNR cuts: snr > $(config["snr_threshold"]), isnr > $(config["snr_inspiral_threshold"])")
    total_snr_not_invc = 0
    for wf_fam in config["waveform_families"]
        println("Waveform family: $(wf_fam)")
        for pno in config["pn_orders"]
            println("  PN order $(pno):")
            println("      Selected events = $(sum(pn_indices[wf_fam][pno]))")
            count_not_invc = sum(noninvertible_after_snr[wf_fam][pno])
            total_snr_not_invc += count_not_invc
            println("      SNR-selected but non-invertible = $(count_not_invc)")
        end

        println("  Final events used across all PN orders: $(sum(summary_indices[wf_fam]))")
        println("  Total SNR-selected but non-invertible across PN orders: $(total_snr_not_invc)")
    end

    shared_index = nothing
    for wf_fam in config["waveform_families"]
        if isnothing(shared_index)
            shared_index = copy(summary_indices[wf_fam])
        else
            shared_index .&= summary_indices[wf_fam]
        end
    end
    println("Events used in all waveform families: $(sum(shared_index))")

    return shared_index
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

Returns
-------
- `Vector{Float64}`: One bootstrap population constraint per realization.
"""
function bootstrapPopConstraint(
    dphi_k::Vector{Float64}, 
    delta_k::Vector{Float64}, 
    n_catalog::Int, 
    n_samples::Int, 
    method::String
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

    constraints = Vector{Float64}(undef, n_samples)
    n_events = length(dphi_k)

    for idx_sample in 1:n_samples
        selected_indices = randperm(n_events)[1:n_catalog]
        constraints[idx_sample] = constraint(
            dphi_k[selected_indices],
            delta_k[selected_indices],
        )
    end

    return constraints
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#
function run_population_analysis(config::Dict)
    fisher_results_file = get_fisher_results_file(config)
    population_results_file = get_population_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    println("Using Fisher results file: $(fisher_results_file)")

    pn_indices = Dict{String, Dict{String, BitVector}}()
    summary_indices = Dict{String, BitVector}()
    noninvertible_after_snr = Dict{String, Dict{String, BitVector}}()
    selected_delta_k = Dict{String, Dict{String, Vector{Float64}}}()
    selected_dphi_k = Dict{String, Dict{String, Vector{Float64}}}()
    phi90 = Dict{String, Dict{String, Vector{Float64}}}()
    bootstrap_constraints = Dict{String, Dict{String, Vector{Float64}}}()

    h5open(fisher_results_file, "r") do file
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
            selected_delta_k[wf_fam] = Dict{String, Vector{Float64}}()
            selected_dphi_k[wf_fam] = Dict{String, Vector{Float64}}()
            phi90[wf_fam] = Dict{String, Vector{Float64}}()
            bootstrap_constraints[wf_fam] = Dict{String, Vector{Float64}}()
            selection_index = summary_indices[wf_fam]

            for pno in config["pn_orders"]
                pno_grp = createSED.pnoString(pno)
                wf_grp = file[pno_grp][config["network"]][wf_fam]

                delta_k = read(wf_grp, "delta_k")
                dphi_k = read(wf_grp, "dphi_k")

                selected_delta_k[wf_fam][pno] = delta_k[selection_index]
                selected_dphi_k[wf_fam][pno] = dphi_k[selection_index]

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
                    )
            end
        end
    end

    shared_index = print_selection_summary(config, pn_indices, summary_indices, noninvertible_after_snr)
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
    )
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
