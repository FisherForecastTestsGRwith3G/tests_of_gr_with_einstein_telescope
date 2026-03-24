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

include("../create_single_event_datasets/createSED.jl")

const par_results_2021 = Dict(
    "GWTC-3 (SEOB)" => Dict(
        "-1" => 0.0007433108784460107,
        "0" => 0.059613216820424045,
        "0.5" => 0.15691276862448023,
        "1" => 0.10345751156561975,
        "1.5" => 0.0811992172352535,
        "2" => 0.3697143089034979,
        "log(2.5)" => 0.27499346193039875,
        "3" => 0.26217884551217885,
        "log(3.)" => 1.9309196583971364,
        "3.5" => 1.0889741092443792,
    ),
    "GWTC-2 (SEOB)" => Dict(
        "-1" => 0.0016434903371840342,
        "0" => 0.04340030721412105,
        "0.5" => 0.13273543814084343,
        "1" => 0.09217054892730547,
        "1.5" => 0.07012688364039743,
        "2" => 0.35655249844439024,
        "log(2.5)" => 0.21647548449350246,
        "3" => 0.3157414170927684,
        "log(3.)" => 1.6407285664042424,
        "3.5" => 0.9428685442198953,
    ),
    "GWTC-2 (Phenom)" => Dict(
        "-1" => 0.002310102895688485,
        "0" => 0.05092705317930542,
        "0.5" => 0.17163379595812028,
        "1" => 0.11436571706841964,
        "1.5" => 0.06776596416236075,
        "2" => 0.5219746773800826,
        "log(2.5)" => 0.1627235343451559,
        "3" => 0.38625337048760455,
        "log(3.)" => 2.130429979529079,
        "3.5" => 0.7764994724454182,
    ),
    "GW170817 (SEOBNRT)" => Dict(
        "-1" => 1.80508774744172e-05,
        "0" => 0.2979771588780663,
        "0.5" => 0.07282102339581283,
        "1" => 0.1046070424050298,
        "1.5" => 0.24357957783944417,
        "2" => 2.9563380093783644,
        "log(2.5)" => 0.7510916171980936,
        "3" => 1.70824386324803,
        "log(3.)" => 7.1208107820649875,
        "3.5" => 8.647729239713513,
    ),
    "GW170817 (PhenomPNRT)" => Dict(
        "-1" => 1.910056446625911e-05,
        "0" => 0.3421820758027727,
        "0.5" => 0.0798284738694681,
        "1" => 0.11519672569317,
        "1.5" => 0.23998773601361828,
        "2" => 3.1910306819759464,
        "log(2.5)" => 0.877695508873479,
        "3" => 2.153824016616575,
        "log(3.)" => 7.450745150499381,
        "3.5" => 8.797346035298109,
    ),
)

const par_results_2025 = Dict(
    "GWTC-3 (SEOB)" => Dict(
        "-1" => 0.0007433108784460107,
        "0" => 0.059613216820424045,
        "0.5" => 0.15691276862448023,
        "1" => 0.10345751156561975,
        "1.5" => 0.0811992172352535,
        "2" => 0.3697143089034979,
        "log(2.5)" => 0.27499346193039875,
        "3" => 0.26217884551217885,
        "log(3.)" => 1.9309196583971364,
        "3.5" => 1.0889741092443792,
    ),
    "GWTC-2 (SEOB)" => Dict(
        "-1" => 0.0016434903371840342,
        "0" => 0.04340030721412105,
        "0.5" => 0.13273543814084343,
        "1" => 0.09217054892730547,
        "1.5" => 0.07012688364039743,
        "2" => 0.35655249844439024,
        "log(2.5)" => 0.21647548449350246,
        "3" => 0.3157414170927684,
        "log(3.)" => 1.6407285664042424,
        "3.5" => 0.9428685442198953,
    ),
    "GWTC-2 (Phenom)" => Dict(
        "-1" => 0.002310102895688485,
        "0" => 0.05092705317930542,
        "0.5" => 0.17163379595812028,
        "1" => 0.11436571706841964,
        "1.5" => 0.06776596416236075,
        "2" => 0.5219746773800826,
        "log(2.5)" => 0.1627235343451559,
        "3" => 0.38625337048760455,
        "log(3.)" => 2.130429979529079,
        "3.5" => 0.7764994724454182,
    ),
    "GW170817 (SEOBNRT)" => Dict(
        "-1" => 1.80508774744172e-05,
        "0" => 0.2979771588780663,
        "0.5" => 0.07282102339581283,
        "1" => 0.1046070424050298,
        "1.5" => 0.24357957783944417,
        "2" => 2.9563380093783644,
        "log(2.5)" => 0.7510916171980936,
        "3" => 1.70824386324803,
        "log(3.)" => 7.1208107820649875,
        "3.5" => 8.647729239713513,
    ),
    "GW170817 (PhenomPNRT)" => Dict(
        "-1" => 1.910056446625911e-05,
        "0" => 0.3421820758027727,
        "0.5" => 0.0798284738694681,
        "1" => 0.11519672569317,
        "1.5" => 0.23998773601361828,
        "2" => 3.1910306819759464,
        "log(2.5)" => 0.877695508873479,
        "3" => 2.153824016616575,
        "log(3.)" => 7.450745150499381,
        "3.5" => 8.797346035298109,
    ),
)

const RESULT_KEYS = collect(keys(par_results_2021))

function get_pn_orders(result_key::AbstractString)
    result_2021 = get(par_results_2021, result_key, nothing)
    result_2025 = get(par_results_2025, result_key, nothing)
    result_2021 === nothing && throw(ArgumentError("Unknown result: $(result_key)"))
    result_2025 === nothing && throw(ArgumentError("Unknown result: $(result_key)"))

    all_orders = union(collect(keys(result_2021)), collect(keys(result_2025)))
    return filter(pno -> pno in all_orders, createSED.pno_list)
end

function make_comparison_plot(result_key::AbstractString)
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
