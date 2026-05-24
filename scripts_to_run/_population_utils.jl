function get_fisher_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "fisher_results_$(config["catalog_tag"]).h5")
end

function get_population_results_file(config::Dict)
    return joinpath(@__DIR__, config["outdir"], "population_results_$(config["bootstrap_tag"]).h5")
end

function derive_population_seed(base_seed::Integer, wf_idx::Integer, pno_idx::Integer, stream_idx::Integer)
    return Int(base_seed + 10_000 * (wf_idx - 1) + 100 * (pno_idx - 1) + stream_idx)
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
