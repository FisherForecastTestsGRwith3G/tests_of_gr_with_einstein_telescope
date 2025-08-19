using Pkg
using HDF5
import JSON

parentdir = dirname(pwd())
Pkg.activate(string(parentdir)*"/GW.jl")
using GW

# include required scripts 
include("_setup_networks.jl")
include("_setup_gr_deviations.jl")
include("_parse_config.jl")

# check if we should run the simulation. Default is to run it, so that first-time users will have no issues running the code
needToEvaluateFisherSNRs = true
if length(ARGS) > 0
    if ARGS[1] == "1"
        needToEvaluateFisherSNRs = true
        println("Fisher matrices and SNRs will be evaluated!")
    elseif ARGS[1] == "0"
        needToEvaluateFisherSNRs = false
        println("Fisher matrices and SNRs will not be evaluated, and instead will be loaded from disk (saves time if you already evaluated them before)!")
    else
        @warn "There was a first input argument handed. But it was $(ARGS[1]) and not \"1\" or \"0\" and thus ignored"
    end
end

# obtain config file from input
user_configs = getUserConfigs()
config_file_name = "config_files/"*user_configs["default_config"]
if length(ARGS) > 1
    config_file_name = ARGS[2]
else
    println("No config file handed, using default!")
end
println("Using config file: $(config_file_name)\n")

configs, simulation_tag = readConfigForA(config_file_name)

# create folder paths
output_folder_name = user_configs["path_output"]*simulation_tag*"/"

#Read catalog and create GR deviations

println("\nRead catalog and calculate GR-deviations")
gr_parameter = ReadCatalog(configs["catalog_name"], folder=user_configs["path_catalog"])

n_events = configs["n_events"]
println("Number of events: ", n_events)

pn_deviation = deltaPnNormal(
    gr_parameter[1][1:n_events],
    gr_parameter[2][1:n_events],
    gr_parameter[3][1:n_events],
    gr_parameter[4][1:n_events],
    gr_parameter[5][1:n_events],
    gr_parameter[6][1:n_events],
    gr_parameter[7][1:n_events],
    gr_parameter[8][1:n_events],
    gr_parameter[9][1:n_events],
    gr_parameter[10][1:n_events],
    gr_parameter[11][1:n_events],
    mu=configs["mu"],
    sigma = configs["sigma"]
    )

gr_deviation_dict = Dict(
    "-1"       => pn_deviation[1],
    "0"        => pn_deviation[2],
    "0.5"      => pn_deviation[3],
    "1"        => pn_deviation[4],
    "1.5"      => pn_deviation[5],
    "2"        => pn_deviation[6],
    "log(2.5)" => pn_deviation[7],
    "3"        => pn_deviation[8],
    "log(3.)"  => pn_deviation[9],
    "3.5"      => pn_deviation[10],
)

### Set up networks 
println("Loading networks:")
networks = Dict()
for nn in configs["network_list"]
    networks[nn] = getNetwork(nn)
    println(nn)
end
println("Finished loading networks.\n")

### Get done the calculations ###
file_name = output_folder_name*"catalog_w_deviations.h5"

# If the user required to not evaluate the Fisher matrices and SNRs, then I will check that the previously evaluated SNR and Fisher
# were calculated using the same PN deviations as the ones just now generated, otherwise everything below would be inconsistent!
if !needToEvaluateFisherSNRs
    # check if the file exists
    if !isfile(file_name)
        throw(ArgumentError("The file $(file_name) does not exist, but you set needToEvaluateFisherSNRs to false. Please set it to true to perform the first run."))
    end

    # check if the PN deviations are the same as the ones in the file
    h5open(file_name, "r") do file
        param_group = file["parameter"]
        for pno in configs["pn_waveforms"]
            saved_deviation = read(param_group, "pn_" * pn_order_dic[pno][2])
            if !all(saved_deviation[1:n_events] .== gr_deviation_dict[pno][1:n_events])
                throw(ArgumentError("The PN deviations in the file $(file_name) do not match the ones just generated. Please set needToEvaluateFisherSNRs to true to re-run the SNR and Fisher evaluations with the new requested injections of the beyond GR deviations!"))
            end
        end
    end
end

println("Writing parameter to $(file_name)")
mkpath(output_folder_name)
h5open(file_name, "w") do file
    param_group = create_group(file, "parameter")
    write(param_group, "mc", gr_parameter[1][1:n_events])
    write(param_group, "eta", gr_parameter[2][1:n_events])
    write(param_group, "chi1", gr_parameter[3][1:n_events])  
    write(param_group, "chi2", gr_parameter[4][1:n_events])  
    write(param_group, "dL", gr_parameter[5][1:n_events])
    write(param_group, "theta", gr_parameter[6][1:n_events])
    write(param_group, "phi", gr_parameter[7][1:n_events])
    write(param_group, "iota", gr_parameter[8][1:n_events])
    write(param_group, "psi",  gr_parameter[9][1:n_events])
    write(param_group, "phiCoal", gr_parameter[10][1:n_events])
    write(param_group, "tcoal", gr_parameter[11][1:n_events])
    write(param_group, "lambda1", gr_parameter[12][1:n_events])
    write(param_group, "lambda2", gr_parameter[13][1:n_events])

    for pno in configs["pn_waveforms"]
        write(param_group, "pn_"*pn_order_dic[pno][2], gr_deviation_dict[pno])
    end
end
    
### calculate fishers ##########################################################
for nn in keys(networks)

    snr_index_global = ones(Bool, n_events)
    inspiral_snr_index_global = ones(Bool, n_events)
    fisher_index_global = ones(Bool, n_events)

    for pno in configs["pn_waveforms"]
        
        pno_name = "pn_"*pn_order_dic[pno][2]
        println("\nProcessing: "*pno_name)
        
        folder_name = output_folder_name * "data/" * nn * "/" * pno_name *"/"
        println("\nStoring results into:")
        println(folder_name)
        mkpath(folder_name)

        #extract paramter
        mc = gr_parameter[1][1:n_events]
        η = gr_parameter[2][1:n_events]
        χ_1 = gr_parameter[3][1:n_events]
        χ_2 = gr_parameter[4][1:n_events]
        dL = gr_parameter[5][1:n_events]
        θ = gr_parameter[6][1:n_events]
        ϕ = gr_parameter[7][1:n_events]
        iota = gr_parameter[8][1:n_events]
        ψ = gr_parameter[9][1:n_events]
        tcoal = gr_parameter[10][1:n_events]
        Φ_coal = gr_parameter[11][1:n_events]
        delta_pn = gr_deviation_dict[pno]

        #choose waveform
        wf = nothing 
        if configs["use_hm_waveform"]
            if configs["use_spinless_waveform"]
                wf = PhenomHM_TIGER_spinless(pn_order_dic[pno][1])
            else
                wf = PhenomHM_TIGER(pn_order_dic[pno][1])
            end
        else
            if configs["use_spinless_waveform"]
                wf = PhenomD_TIGER_spinless(pn_order_dic[pno][1])
            else
                wf = PhenomD_TIGER(pn_order_dic[pno][1])
            end
        end

        network = networks[nn]

        #calculate fisher
        if needToEvaluateFisherSNRs

            println("Calculating Fisher matrices and SNRs")
            @time fisher_matrices, snrs = FisherMatrix(
                wf,
                network,
                mc, 
                η, 
                χ_1, 
                χ_2, 
                dL, 
                θ, 
                ϕ, 
                iota, 
                ψ, 
                tcoal, 
                Φ_coal, 
                delta_pn, 
                auto_save=false, 
                return_SNR=true, 
                useEarthMotion=true
            )

            println("Calculating inspiral SNRs")
            #I evaluate the SNR at the end of the "inspiral" phase, as defined in the Phenom waveform models
            f_inspiral_cutoff = @. 0.018 / (  mc / η^(3. /5.) ) / GMsun_over_c3

            @time inspiral_snrs = SNR(
                wf,
                network,
                mc, 
                η, 
                χ_1, 
                χ_2, 
                dL, 
                θ, 
                ϕ, 
                iota, 
                ψ, 
                tcoal, 
                delta_pn, 
                auto_save=false, 
                useEarthMotion=true,
                fmax = f_inspiral_cutoff
            )
            
        else
            
            
            println("\nSkipping the simulation since data should already exist")
            println("Loading the results from outputfolder: ")
            print(folder_name)
            print()

            # read fisher matrices (which are already stored)
            filename = folder_name * "fishers.h5"
            fisher_matrices = h5open(filename, "r") do file
                read(file, "matrices")  
            end
    
            # read snrs (which are already stored)
            filename = folder_name * "snrs.h5"
            snrs = h5open(filename, "r") do file
                read(file, "values")  
            end
            
            # read inspiral snrs (which are already stored)
            filename = folder_name * "inspiral_snrs.h5"
            inspiral_snrs = h5open(filename, "r") do file
                read(file, "values")  
            end

            # Select only the first n_events, in case the fisher matrices and SNRs were calculated for more events than requested
            fisher_matrices = fisher_matrices[1:n_events, :, :]
            snrs = snrs[1:n_events]
            inspiral_snrs = inspiral_snrs[1:n_events]
            
        end

        ### postprocessing #########################################################

        # invert fisher matrices and calculate errors 
        covariance_matrices = CovMatrix(fisher_matrices)

        deltak_pn = Array{Float64}(undef, n_events)
        fisher_inverted = Array{Bool}(undef, n_events)
        not_inverted = 0

        for idx_event = 1:n_events
            cmat = covariance_matrices[idx_event, :,:]
            if all(cmat .== 0.0)
                # TODO: Find a better way to deal with non-invertible fisher matrices
                deltak_pn[idx_event] = NaN
                fisher_inverted[idx_event] = false
                not_inverted += 1
            else
                # take the squareroot of the 12th parameter. This is the GR deviation.
                deltak_pn[idx_event] = sqrt(cmat[12, 12])
                fisher_inverted[idx_event] = true
            end
        end

        print("$(not_inverted) Fisher matrices (out of $(n_events)) could not be inverted, or had a SNR below the threshold.")

        # caclulate the expected deviations 
        # TODO: Do a better calculation than this. This is provisorical
        #       See also: https://arxiv.org/abs/gr-qc/0703086 equation (35) for instructions
        dphi0k_pn = gr_deviation_dict[pno] .+ deltak_pn .* randn(n_events)

        # calculate indices and update global indices
        snr_index = convert(Vector{Bool}, snrs .> configs["snr_thresh"])
        inspiral_snr_index = convert(Vector{Bool}, inspiral_snrs .> configs["inspiral_snr_thresh"])
        snr_index_global = convert(Vector{Bool}, snr_index .& snr_index_global)
        inspiral_snr_index_global = convert(Vector{Bool}, inspiral_snr_index .& inspiral_snr_index_global)
        fisher_index_global = convert(Vector{Bool}, fisher_inverted .& fisher_index_global)

        # save everything to .h5 file 
        filename = folder_name * "single_event_measurements.h5"
        h5open(filename, "w") do file
            write(file, "dphit_k", gr_deviation_dict[pno])
            write(file, "dphi0_k", dphi0k_pn)
            write(file, "delta_k", deltak_pn)
        end

        filename = folder_name * "fishers.h5"
        h5open(filename, (isfile(filename) ? "r+" : "w")) do file
            if needToEvaluateFisherSNRs
                if haskey(file, "matrices")
                    HDF5.delete_object(file, "matrices")
                end
                write(file, "matrices", fisher_matrices)
            end
            if haskey(file, "index")
                HDF5.delete_object(file, "index")
            end
            write(file, "index", fisher_inverted)
        end

        filename = folder_name * "snrs.h5"
        h5open(filename, (isfile(filename) ? "r+" : "w")) do file
            if needToEvaluateFisherSNRs
                if haskey(file, "values")
                    HDF5.delete_object(file, "values")
                end
                write(file, "values", snrs)
            end
            if haskey(file, "index")
                HDF5.delete_object(file, "index")
            end
            write(file, "index", snr_index)
        end

        filename = folder_name * "inspiral_snrs.h5"
        h5open(filename, (isfile(filename) ? "r+" : "w")) do file
            if needToEvaluateFisherSNRs
                if haskey(file, "values")
                    HDF5.delete_object(file, "values")
                end
                write(file, "values", inspiral_snrs)
            end
            if haskey(file, "index")
                HDF5.delete_object(file, "index")
            end
            write(file, "index", inspiral_snr_index)
        end

    end

    # saving global indices
    
    # For additional flexibility in the following scripts, I save both total_index_global-like variables, assuming both use_inspiral_snr_thresh=true and use_inspiral_snr_thresh=false anyway
    total_index_global_with_inspiral_snr_cut = convert(Vector{Bool}, snr_index_global .& fisher_index_global .& inspiral_snr_index_global)
    total_index_global_without_inspiral_snr_cut = convert(Vector{Bool}, snr_index_global .& fisher_index_global)

    # But in the end, I set the total_global_index variables, which is the one automatically used in the other scripts, according to the current setting of use_inspiral_snr_thresh from the config file
    if configs["use_inspiral_snr_thresh"]
        total_index_global = total_index_global_with_inspiral_snr_cut
    else
        total_index_global = total_index_global_without_inspiral_snr_cut
    end

    println("\nStatistics for network $(nn):")
    println("Total number of events: ", length(total_index_global))
    println("Total number of events after SNR cut (snrs > ", configs["snr_thresh"] , ")): ", sum(snr_index_global), " (", round(sum(snr_index_global) * 100. / length(total_index_global) , digits=2), "%)")
    println("Total number of events after inspiral SNR cut (inspiral_snrs > ", configs["inspiral_snr_thresh"] , "): ", sum(inspiral_snr_index_global), " (", round(sum(inspiral_snr_index_global) * 100. / length(total_index_global), digits=2), "%)")
    println("Total number of events after Fisher cut: ", sum(fisher_index_global), " (", round(sum(fisher_index_global) * 100. / length(total_index_global), digits=2), "%)")
    println("Total number of events after Fisher and SNR cut: ", sum(total_index_global_without_inspiral_snr_cut), " (", round(sum(total_index_global_without_inspiral_snr_cut) * 100. / length(total_index_global) , digits=2), "%)")
    println(",- Total number of events after Fisher and SNR and inspiral SNR cut: ", sum(total_index_global_with_inspiral_snr_cut), " (", round(sum(total_index_global_with_inspiral_snr_cut) * 100. / length(total_index_global), digits=2), "%)")
    if configs["use_inspiral_snr_thresh"]
        println("`-> currently applying inspiral SNR cut")
    else
        println("`-> currently not applying inspiral SNR cut")
    end

    filename = output_folder_name * "data/" * nn * "/global_indices.h5"
    h5open(filename, "w") do file
        write(file, "snr", snr_index_global)
        write(file, "inspiral_snr", inspiral_snr_index_global)
        write(file, "fisher", fisher_index_global)
        write(file, "total", total_index_global)
        write(file, "total_with_inspiral_snr_cut", total_index_global_with_inspiral_snr_cut)
        write(file, "total_without_inspiral_snr_cut", total_index_global_without_inspiral_snr_cut)
    end
end

