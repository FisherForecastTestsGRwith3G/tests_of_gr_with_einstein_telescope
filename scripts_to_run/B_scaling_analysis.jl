using TOML
using HDF5
using Printf
using Random
using Dates
using Statistics

include("_config_parser.jl")
include("_hdf5_metadata.jl")
include("_population_utils.jl")
include("../create_single_event_datasets/createSED.jl")
include("../hierachical_combination/hierDist.jl")

#----------------------------------------------------------------------------#
# Auxiliary functions
#----------------------------------------------------------------------------#

function get_scaling_results_file(config::Dict)
    return joinpath(@__DIR__, config["scaling_outdir"], "scaling_results_$(config["scaling_tag"]).h5")
end

"""
Write scaling-analysis results to an HDF5 file.

The file is created at `scaling_results_file`. Existing files at that path are
overwritten. Results are stored under
`/<network>/<waveform_family>/<pn_order>/`, with one dataset per summary curve.

Arguments
---------
- `scaling_results_file`: Output HDF5 file path.
- `config`              : Analysis configuration dictionary.
- `summary_indices`     : Per-waveform-family event masks used for selection.
- `selected_delta_k`    : Selected `delta_k` values used in the scaling analysis.
- `selected_dphi_k`     : Selected `dphi_k` values used in the scaling analysis.
- `scaling_results`     : Per-waveform-family, per-PN-order scaling summaries.
"""
function write_scaling_results_hdf5(
    scaling_results_file::AbstractString,
    config::Dict,
    summary_indices::Dict{String, BitVector},
    selected_delta_k::Dict{String, Dict{String, Vector{Float64}}},
    selected_dphi_k::Dict{String, Dict{String, Vector{Float64}}},
    scaling_results::Dict{String, Dict{String, Any}},
)
    mkpath(dirname(scaling_results_file))

    h5open(scaling_results_file, "w") do file
        write_top_level_output_metadata!(file)
        attrs(file)["scaling_tag"] = config["scaling_tag"]
        attrs(file)["n_max"]       = config["scaling_n_max"]
        attrs(file)["n_steps"]     = config["scaling_n_steps"]
        attrs(file)["n_sample"]    = config["scaling_n_sample"]
        attrs(file)["method"]      = config["scaling_method"]

        network_group = create_group(file, config["network"])

        for wf_fam in config["waveform_families"]
            wf_grp = create_group(network_group, wf_fam)
            attrs(wf_grp)["n_selected"] = sum(summary_indices[wf_fam])
            write(wf_grp, "summary_indices", Vector{Bool}(summary_indices[wf_fam]))

            for pno in config["pn_orders"]
                result = scaling_results[wf_fam][pno]
                pno_grp = create_group(wf_grp, createSED.pnoString(pno))

                attrs(pno_grp)["n_selected"] = length(selected_delta_k[wf_fam][pno])
                write(pno_grp, "population_sizes", result.population_sizes)

                pop_grp = create_group(pno_grp, "population")
                write(pop_grp, "median", result.population.median)
                write(pop_grp, "q05", result.population.q05)
                write(pop_grp, "q95", result.population.q95)

                best_grp = create_group(pno_grp, "best_single")
                write(best_grp, "median", result.best_single.median)
                write(best_grp, "q05", result.best_single.q05)
                write(best_grp, "q95", result.best_single.q95)
            end
        end
    end
end

"""
Estimate how population constraints scale with the number of events.

For each population size between 10 and `n_max`, the function draws `n_samples`
subsets without replacement from `dphi_k` and `delta_k`, evaluates the requested
population-constraint method on each subset, and records the best single-event
constraint in the same subset. It returns the median, 5%, and 95% quantiles of
both constraint families as functions of population size.

The input catalog must contain enough distinct `n_max`-event combinations that
the probability of sampling the same combination within `n_samples` draws is
below 0.1%.

Arguments
---------
- `dphi_k`   : Mean values of the single-event Gaussian constraints.
- `delta_k`  : Standard deviations of the single-event Gaussian constraints.
- `n_max`    : Largest population size to sample.
- `n_steps`  : Number of population sizes between 10 and `n_max`.
- `n_samples`: Number of random subsets evaluated per population size.
- `seed`     : Random seed used for subset draws.
- `method`   : Population-constraint method. Currently supported: `"prodL"`.

Returns
-------
- Named tuple with `population_sizes`, `population`, and `best_single`, where
  `population` and `best_single` each contain `median`, `q05`, and `q95` vectors.
"""
function scaling_analysis(
    dphi_k::Vector{Float64}, 
    delta_k::Vector{Float64},
    n_max::Int,
    n_steps::Int,
    n_samples::Int = 1000,
    seed::Int=1234,
    method::String="prodL"
)

    # Input parsing checks
    length(dphi_k) == length(delta_k) ||
        throw(ArgumentError("`dphi_k` and `delta_k` must have the same length."))
    n_events = length(dphi_k)
    n_max >= 10 || throw(ArgumentError("`n_max` must be at least 10."))
    n_max <= n_events || throw(ArgumentError("`n_max` cannot exceed the number of available events."))
    n_steps > 0 || throw(ArgumentError("`n_steps` must be positive."))
    n_samples > 0 || throw(ArgumentError("`n_samples` must be positive."))

    n_combinations = binomial(big(n_events), big(n_max))
    min_combinations = BigFloat(n_samples) * BigFloat(n_samples - 1) /
                       (-2 * log1p(BigFloat(-0.001)))
    if !(n_combinations >= min_combinations) 
        throw(ArgumentError(
            "Too few distinct event combinations for scaling analysis: choose($(n_events), $(n_max)) = $(n_combinations), giving collision probability above 0.1% for $(n_samples) draws."
        ))
    end

    # select method
    if method == "prodL"
        constraint(mu, sig) = prodLPopDeviationConstraint(mu, sig)
    else
        throw(ArgumentError("Unknown bootstrap method: $(method)"))
    end

    # perform scaling analysis
    population_sizes = round.(Int, range(10, n_max; length=n_steps))
    population_constraints = Matrix{Float64}(undef, n_steps, n_samples)
    best_single_constraints = Matrix{Float64}(undef, n_steps, n_samples)
    single_event_constraints = HierDist.estimate0Symmetric90CiGaussian(dphi_k, delta_k)

    rng = MersenneTwister(seed)
    permutation = collect(1:n_events)

    for (idx_pop, k_pop) in pairs(population_sizes)
        selected_dphi = Vector{Float64}(undef, k_pop)
        selected_delta = Vector{Float64}(undef, k_pop)

        for idx_sample in 1:n_samples
            best_constraint = Inf

            # Samples k_pop distinct event indices without needing to shuffle 
            # all n_events using a Fisher-Yates shuffle
            @inbounds for idx_event in 1:k_pop
                swap_index = rand(rng, idx_event:n_events)
                permutation[idx_event], permutation[swap_index] =
                    permutation[swap_index], permutation[idx_event]
                selected_index = permutation[idx_event]
                selected_dphi[idx_event] = dphi_k[selected_index]
                selected_delta[idx_event] = delta_k[selected_index]
                event_constraint = single_event_constraints[selected_index]
                if event_constraint < best_constraint
                    best_constraint = event_constraint
                end
            end

            population_constraints[idx_pop, idx_sample] = constraint(selected_dphi, selected_delta)
            best_single_constraints[idx_pop, idx_sample] = best_constraint
        end
    end

    pop_median = Vector{Float64}(undef, n_steps)
    pop_q05 = Vector{Float64}(undef, n_steps)
    pop_q95 = Vector{Float64}(undef, n_steps)
    best_median = Vector{Float64}(undef, n_steps)
    best_q05 = Vector{Float64}(undef, n_steps)
    best_q95 = Vector{Float64}(undef, n_steps)

    for idx_pop in 1:n_steps
        pop_q05[idx_pop], pop_median[idx_pop], pop_q95[idx_pop] =
            quantile(@view(population_constraints[idx_pop, :]), [0.05, 0.50, 0.95])
        best_q05[idx_pop], best_median[idx_pop], best_q95[idx_pop] =
            quantile(@view(best_single_constraints[idx_pop, :]), [0.05, 0.50, 0.95])
    end

    return (
        population_sizes=population_sizes,
        population=(median=pop_median, q05=pop_q05, q95=pop_q95),
        best_single=(median=best_median, q05=best_q05, q95=best_q95),
    )
end

#----------------------------------------------------------------------------#
# MAIN FUNCTION
#----------------------------------------------------------------------------#

"""
Run the scaling analysis and write its results to disk.

The function loads Fisher-analysis outputs, builds the shared selection masks
with `collect_population_inputs`, applies each waveform family's
`summary_indices` mask to every configured PN order, and evaluates
`scaling_analysis` on those selected samples.
"""
function run_scaling_analysis(config::Dict; verbose::Bool=true)
    fisher_results_file  = get_fisher_results_file(config)
    scaling_results_file = get_scaling_results_file(config)
    isfile(fisher_results_file) || throw(ArgumentError("Fisher results file does not exist: $(fisher_results_file)"))

    println("Using Fisher results file: $(fisher_results_file)")

    h5open(fisher_results_file, "r") do file
        _, summary_indices, noninvertible_after_snr, delta_k_data, dphi_k_data =
            collect_population_inputs(file, config)

        selected_delta_k = Dict{String, Dict{String, Vector{Float64}}}()
        selected_dphi_k = Dict{String, Dict{String, Vector{Float64}}}()
        scaling_results = Dict{String, Dict{String, Any}}()

        if verbose
            println("\nScaling analysis")
            println("Applied SNR cuts: snr > $(config["snr_threshold"]), isnr > $(config["snr_inspiral_threshold"])")
            println("Scaling settings: n_max = $(config["scaling_n_max"]), n_steps = $(config["scaling_n_steps"]), n_sample = $(config["scaling_n_sample"])")
        end

        for (idx_wf, wf_fam) in enumerate(config["waveform_families"])
            summary_index = summary_indices[wf_fam]
            selected_delta_k[wf_fam] = Dict{String, Vector{Float64}}()
            selected_dphi_k[wf_fam] = Dict{String, Vector{Float64}}()
            scaling_results[wf_fam] = Dict{String, Any}()

            for (idx_pno, pno) in enumerate(config["pn_orders"])
                selected_delta_k[wf_fam][pno] = delta_k_data[wf_fam][pno][summary_index]
                selected_dphi_k[wf_fam][pno] = dphi_k_data[wf_fam][pno][summary_index]

                scaling_results[wf_fam][pno] = scaling_analysis(
                    selected_dphi_k[wf_fam][pno],
                    selected_delta_k[wf_fam][pno],
                    config["scaling_n_max"],
                    config["scaling_n_steps"],
                    config["scaling_n_sample"],
                    derive_population_seed(config["seed"], idx_wf, idx_pno, 3),
                    config["scaling_method"],
                )

                if verbose
                    @printf(
                        "  %-10s PN %-6s selected: %d/%d, non-invertible after SNR: %d\n",
                        wf_fam,
                        pno,
                        length(selected_delta_k[wf_fam][pno]),
                        length(delta_k_data[wf_fam][pno]),
                        sum(noninvertible_after_snr[wf_fam][pno])
                    )
                end
            end
        end

        write_scaling_results_hdf5(
            scaling_results_file,
            config,
            summary_indices,
            selected_delta_k,
            selected_dphi_k,
            scaling_results,
        )
    end

    println("Scaling results written to $(scaling_results_file)")
    return scaling_results_file
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

    return run_scaling_analysis(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
