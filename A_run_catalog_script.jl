using Pkg
using HDF5

parentdir = dirname(pwd())
Pkg.activate(string(parentdir)*"/GW.jl")
using GW

# include required scripts 
include("_setup_networks.jl")
include("_setup_gr_deviations.jl")
include(".paths.jl")

################################################################################
## Specify simulation specs in this part of the script
# ||||||||
# vvvvvvvv

PhD = "Joachim"
path_catalog, path_output = whoIsThere(PhD)
needToRun = true # set to true if you want to run the simulation
HM = false # set to true if you want to run the simulation with the HM waveform

# how is this simulation called
simulation_tag = "test_mcmc"

# network specs
network_names = ["ETS"]
                #["ETS", "network_0_15km", "network_45_15km"]

# specs of catalog
n_events      = 100
source_type   = "BBH"
catalog_name  = "BGR_TIGER_10k.h5"
pn_orders = ["0", "0.5", "1", "1.5"] 
            #["-1", "0", "0.5", "1", "1.5", "2", "log(2.5)", "3", "log(3.)", "3.5"]   

# snr threshold
snr_thresh = 12.

# specify GR deviations
mu = 0.0
sigma = 0.0025

# ^^^^^^^^
# ||||||||
## Specify simulation specs in this part of the script
################################################################################

output_folder_name = path_output*"output/"*simulation_tag*"/"

pn_order_dic = Dict(
    "-1"       => (-1.0    ,"minus_one"),
    "0"        => (0.0     ,"zero"), 
    "0.5"      => (0.5     ,"half"),
    "1"        => (1.0     ,"one"), 
    "1.5"      => (1.5     ,"one_half"),
    "2"        => (2.0     ,"two"),
    "log(2.5)" => (log(2.5),"log_two_half"),
    "3"        => (3.0     ,"three"), 
    "log(3.)"  => (log(3.) ,"log_three"),
    "3.5"      => (3.5     ,"three_half")
);

### Read a catalog and create GR deviations
println("Read catalog and calculate GR-deviations")
gr_parameter = ReadCatalog(catalog_name, folder=path_catalog)
# pn_deviation = farrEtAl(
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
    mu=mu*ones(n_events),
    sigma = sigma*ones(n_events)
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
println("Collecting networks:")
networks = Dict()
for nn in network_names
    networks[nn] = getNetwork(nn)
    println(nn)
end

### Get done the calculations ###
file_name = output_folder_name*"catalog_w_deviations.h5"

println("Writing parameter to $(file_name)")
mkpath(output_folder_name)
h5open(file_name, "w") do file
    gr_param_group = create_group(file, "gr_parameter")
    write(gr_param_group, "mc", gr_parameter[1][1:n_events])
    write(gr_param_group, "eta", gr_parameter[2][1:n_events])
    write(gr_param_group, "chi1", gr_parameter[3][1:n_events])  
    write(gr_param_group, "chi2", gr_parameter[4][1:n_events])  
    write(gr_param_group, "dL", gr_parameter[5][1:n_events])
    write(gr_param_group, "theta", gr_parameter[6][1:n_events])
    write(gr_param_group, "phi", gr_parameter[7][1:n_events])
    write(gr_param_group, "iota", gr_parameter[8][1:n_events])
    write(gr_param_group, "psi",  gr_parameter[9][1:n_events])
    write(gr_param_group, "phiCoal", gr_parameter[10][1:n_events])
    write(gr_param_group, "tcoal", gr_parameter[11][1:n_events])
    write(gr_param_group, "lambda1", gr_parameter[12][1:n_events])
    write(gr_param_group, "lambda2", gr_parameter[13][1:n_events])

    gr_deviation_group = create_group(file, "gr_deviation")
    for pno in pn_orders
        write(gr_deviation_group, "pn_"*pn_order_dic[pno][2], gr_deviation_dict[pno])
    end
end
    
### calculate fishers ##########################################################
for nn in keys(networks)

    snr_index_global = ones(Bool, n_events)
    fisher_index_global = ones(Bool, n_events)

    for pno in pn_orders
        
        pno_name = "pn_"*pn_order_dic[pno][2]
        println("\n"*"#"^81)
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

        #calculate fisher
        wf = nothing 

        if HM
            wf = PhenomHM_TIGER(pn_order_dic[pno][1])
        else
            wf = PhenomD_TIGER(pn_order_dic[pno][1])
        end

        network = networks[nn]

        if needToRun
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
        else #Load the matrices instead of rerunning the analysis

            
            println("\nSkipping the simulation since data should already exist")
            println("Loading the results from outputfolder: ")
            print(folder_name)
            print()

            # read fisher matrices (which are already stored)
            filename = folder_name * "fishers.h5"
            fisher_matrices = h5open(filename, "r") do file
                fisher_matrices = read(file, "matrices")  
            end
    
            # read snrs (which are already stored)
            filename = folder_name * "snrs.h5"
            snrs = h5open(filename, "r") do file
                snrs = read(file, "values")  
            end

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

        print("$(not_inverted) Fisher matrices could not be inverted")

        # caclulate the expected deviations 
        # TODO: Do a better calculation than this. This is provisorical
        #       See also: https://arxiv.org/abs/gr-qc/0703086 equation (35) for instructions
        dphi0k_pn = gr_deviation_dict[pno] .+ deltak_pn .* randn(n_events)

        # calculate indices and update global indices
        snr_index = convert(Vector{Bool}, snrs .> snr_thresh)
        snr_index_global = convert(Vector{Bool}, snr_index .& snr_index_global)
        fisher_index_global = convert(Vector{Bool}, fisher_inverted .& fisher_index_global)

        # save everything to .h5 file 
        filename = folder_name * "single_event_measurements.h5"
        h5open(filename, "w") do file
            write(file, "dphit_k", gr_deviation_dict[pno])
            write(file, "dphi0_k", dphi0k_pn)
            write(file, "delta_k", deltak_pn)
        end

        filename = folder_name * "fishers.h5"
        h5open(filename, "w") do file
            write(file, "matrices", fisher_matrices)
            write(file, "index", fisher_inverted)
        end

        filename = folder_name * "snrs.h5"
        h5open(filename, "w") do file
            write(file, "values", snrs)
            write(file, "index", snr_index)
        end
    end

    # saving global index
    total_index_global = convert(Vector{Bool}, snr_index_global .& fisher_index_global)

    filename = output_folder_name * "data/" * nn * "/global_indices.h5"
    h5open(filename, "w") do file
        write(file, "snr", snr_index_global)
        write(file, "fisher", fisher_index_global)
        write(file, "total", total_index_global)
    end
end


