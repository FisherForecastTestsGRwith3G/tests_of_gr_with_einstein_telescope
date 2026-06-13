using HDF5
using Trapz
using KernelDensity, Statistics
using Serialization
using Base.Threads
using Interpolations

# include required scripts
include("../hierachical_combination/hierarchical_distribution.jl")
include("../hierachical_combination/stat_utils.jl")
include("_config_parser_grid.jl")
include("_grid_utils.jl")
include("../scripts_to_run/_config_parser.jl")
include("../create_single_event_datasets/createSED.jl")
include("../scripts_to_run/_population_utils.jl")
using .createSED: pnoString

const N_POINTS_HYPER = 2000
const K_SPREAD_HYPER = 30.0

function wrapper_3sigma_local(n_events_used, dphi0_k, delta_k, center_mu, center_sig)
    n_events_used = Int(round(n_events_used))
    dphi0_k = dphi0_k[1:n_events_used]
    delta_k = delta_k[1:n_events_used]
    
    number_of_zeros = count(x -> x == 0.0, delta_k)
    if number_of_zeros > 0
        @warn "There are $(number_of_zeros) zero values in delta_k. This may lead to issues in the hyperparameter distribution calculation."
    end


    spread = sqrt(1.0 / sum(1.0 ./ delta_k.^2))

    mu_limit = (center_mu - K_SPREAD_HYPER * spread, center_mu + K_SPREAD_HYPER * spread)
    sig_limit = (max(0.0, center_sig - K_SPREAD_HYPER * spread), center_sig + K_SPREAD_HYPER * spread)

    mu_values = collect(LinRange(mu_limit[1], mu_limit[2], N_POINTS_HYPER))
    sig_values = collect(LinRange(sig_limit[1], sig_limit[2], N_POINTS_HYPER))

    p_mu_sig, _, _, _, _ = getDistributionOnGrid(mu_values, sig_values, dphi0_k, delta_k)

    itp = interpolate((mu_values, sig_values), p_mu_sig, Gridded(Linear()))
    p_mu_sig_interp = extrapolate(itp, 0.0)

    local level, p_GR
    try
        level = getContourLevelOGD(p_mu_sig, 0.9889)
        p_GR = p_mu_sig_interp(0.0, 0.0)
    catch
        return 1.0
    end

    # Positive means the GR point lies outside the 3-sigma contour.
    return level - p_GR
end

# obtain config file from input
config_file_name = ARGS[1]
println("Using config file: $(config_file_name)\n")

configs = read_config_grid(config_file_name)

mu_vec = configs["mu_vec"]
sigma_vec = configs["sigma_vec"]
pn_string = configs["pn_orders"][1]
pn_tag = pnoString(pn_string)
network = configs["network"]
waveform = configs["waveform_families"][1]
run_tag = grid_run_tag(network, waveform; grid_tag=configs["grid_tag"])
n_median = configs["n_median"]
grid_size = length(mu_vec) * length(sigma_vec)

output_folder_name = abspath(joinpath(@__DIR__, configs["outdir"], "grid", network, pn_tag))
plot_folder_name = abspath(joinpath(@__DIR__, configs["plot_outdir"], "grid", network, pn_tag))

mkpath(output_folder_name)

println("\nAnalyzing grid composed of $(grid_size) grid points\n")
println("\nCalculating $(n_median) median values for each grid point\n")

idx = 0
res = fill(NaN, length(mu_vec), length(sigma_vec))
for i in eachindex(mu_vec)
    for j in eachindex(sigma_vec)
        global idx += 1
        println("\nAnalyzing data for index $(idx) out of $(grid_size)")

        if mu_vec[i] == 0.0 && sigma_vec[j] == 0.0
            continue
        end

        mu = mu_vec[i]
        sigma = sigma_vec[j]
        header = grid_point_header(pn_tag, idx, mu, sigma)

        println("Mu: ", mu, " Sigma: ", sigma)

        grid_config_file = joinpath(@__DIR__, "config_files", "grid", run_tag, "config_A_$(header).toml")
        if !isfile(grid_config_file)
            @warn "Skipping missing grid config: $(grid_config_file)"
            continue
        end
        grid_config = read_config(grid_config_file)
        fisher_file = get_fisher_results_file(grid_config)

        if !isfile(fisher_file)
            @warn "Skipping missing Fisher file: $(fisher_file)"
            continue
        end

        dphi0_k, delta_k = h5open(fisher_file, "r") do file
            if !haskey(file, pn_tag)
                @warn "Skipping incomplete Fisher file without PN group $(pn_tag): $(fisher_file)"
                return Float64[], Float64[]
            end
            _, summary_indices, _, delta_k_data, dphi_k_data = collect_population_inputs(file, grid_config)
            return dphi_k_data[waveform][pn_string][summary_indices[waveform]],
                   delta_k_data[waveform][pn_string][summary_indices[waveform]]
        end
        isempty(dphi0_k) && continue
        println("Number of valid events for this grid point: ", length(dphi0_k))

        res[i, j] = reshuffling_bisection(n_median, wrapper_3sigma_local, dphi0_k, delta_k, mu, sigma)[1]
    end
end

h5open(joinpath(output_folder_name, "results.h5"), "w") do file
    write(file, "results", res)
    write(file, "mu_vec", mu_vec)
    write(file, "sigma_vec", sigma_vec)
    write(file, "pn_string", pn_string)
end
println("\nSaved intermediate grid data in: $(joinpath(output_folder_name, "results.h5"))")
println("Run C_grid_plot.jl with the same config to generate the figure.")
