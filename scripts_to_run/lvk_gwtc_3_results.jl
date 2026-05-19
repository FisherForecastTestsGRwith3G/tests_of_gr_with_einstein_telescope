# The LVK paper: Tests of General Relativity with GWTC-3 (https://arxiv.org/abs/2112.06861) 
# was updated in 2025. In this script we compare the results of the two different releases from  
# the publicly released data.
#
# Original data release: 2021
#   https://dcc.ligo.org/LIGO-P2100456-v1
#
# Updated data release: 2025
#   https://dcc.ligo.org/LIGO-P2100456-v2
#

using Plots
using LaTeXStrings
import JSON

include("../create_single_event_datasets/createSED.jl")

const RESULTS_2021_JSON = joinpath(@__DIR__, "lvk_gwtc_3_results_2021.json")
const RESULTS_2025_JSON = joinpath(@__DIR__, "lvk_gwtc_3_results_2025.json")

function load_par_results(json_path::AbstractString)
    return JSON.parsefile(json_path; dicttype=Dict{String, Any})
end

function load_all_par_results()
    return load_par_results(RESULTS_2021_JSON), load_par_results(RESULTS_2025_JSON)
end

const RESULT_KEYS = collect(keys(load_par_results(RESULTS_2021_JSON)))

function get_pn_orders(result_key::AbstractString)
    par_results_2021, par_results_2025 = load_all_par_results()

    result_2021 = get(par_results_2021, result_key, nothing)
    result_2025 = get(par_results_2025, result_key, nothing)
    result_2021 === nothing && throw(ArgumentError("Unknown result: $(result_key)"))
    result_2025 === nothing && throw(ArgumentError("Unknown result: $(result_key)"))

    all_orders = union(collect(keys(result_2021)), collect(keys(result_2025)))
    return filter(pno -> pno in all_orders, createSED.pno_list)
end

function make_comparison_plot(result_key::AbstractString)
    par_results_2021, par_results_2025 = load_all_par_results()
    pn_orders = get_pn_orders(result_key)
    values_2021 = [par_results_2021[result_key][pno] for pno in pn_orders]
    values_2025 = [par_results_2025[result_key][pno] for pno in pn_orders]
    xticks = collect(1:length(pn_orders))
    labels = createSED.pnoLatex.(pn_orders)

    plt = plot(
        xticks,
        values_2021;
        label="2021 release",
        lw=2.5,
        marker=:circle,
        markersize=7,
        color="#1f77b4",
        xlabel="PN order",
        ylabel=L"\delta \phi_k",
        yscale=:log10,
        xticks=(xticks, labels),
        title=result_key,
        legend=:topright,
        framestyle=:box,
        grid=true,
        minorgrid=true,
        size=(1100, 700),
        dpi=200,
    )

    plot!(
        plt,
        xticks,
        values_2025;
        label="2025 release",
        lw=2.5,
        marker=:diamond,
        markersize=7,
        color="#d62728",
    )

    return plt
end

function main(args=ARGS)
    par_results_2021, par_results_2025 = load_all_par_results()

    if length(args) != 1
        script_name = basename(@__FILE__)
        valid_results = join(sort(RESULT_KEYS), ", ")
        throw(ArgumentError("Usage: julia $(script_name) \"<result>\"\nValid results: $(valid_results)"))
    end

    result_key = args[1]
    haskey(par_results_2021, result_key) || throw(ArgumentError("Unknown result: $(result_key)"))
    haskey(par_results_2025, result_key) || throw(ArgumentError("Missing 2025 data for result: $(result_key)"))

    plt = make_comparison_plot(result_key)
    display(plt)
    gui(plt)
    println("Plot window opened for $(result_key). Press Enter to close the script.")
    try
        readline()
    catch
    end

    return plt
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
