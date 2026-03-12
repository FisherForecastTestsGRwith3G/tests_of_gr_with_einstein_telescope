using Pkg
using HDF5
import JSON

project_dir = dirname(pwd())
parentdir = dirname(project_dir)
Pkg.activate(joinpath(parentdir,"GW.jl"))
using GW

include(joinpath(project_dir, "supplementary/_utils.jl"))

include(joinpath(project_dir, "_setup_networks.jl"))

output_folder_name = joinpath(project_dir, "supplementary/data")

## DESCRIPTION
################################################################################
# This script recalculates SNRs and Fischer matrices from a specified catalog 
# for the LVK O3 network, and stores them in a data directory.       
################################################################################
## Specify script parameter here
#
n_events = 100000#000
fmin     = 10
fmax     = 4000
#
##
################################################################################

# load network
network = getNetwork("LHV_O3")

# load parameter
parameters = ReadCatalog(
    "BGR_TIGER_1M.h5",
    folder = joinpath(project_dir, "catalogs/"),
    redshift = true
    )

# deploy a redshift cutoff
redshift = parameters[14][:]
idx_redshift = redshift .< 1 #Speeds things up a lot and most events are below this threshold anyway!

# deploy the n_events cutoff
n_available = sum(idx_redshift)
if n_available > n_events
    idx_redshift = idx_redshift .& (cumsum(idx_redshift) .<= n_events)
else
    n_events = n_available
end 

# select parameter and organize them as matrix
parameters_mat = zeros(14, n_events)
for i in 1:14
    parameters_mat[i, 1:n_events] = parameters[i][idx_redshift]
end

# set up waveforms
waveforms = Dict(
    #"PhenomD"   => PhenomD(),
    #"PhenomHM"  => PhenomHM(),
    "PhenomD_pn_log_minus_one" => PhenomD_TIGER_spinless(-1.),
    #"PhenomD_pn_log_two_half"  => PhenomD_TIGER_spinless(log(2.5)),
    #"PhenomD_pn_three"         => PhenomD_TIGER_spinless(3.),
    "PhenomD_pn_three_half"    => PhenomD_TIGER_spinless(3.5),
    "PhenomHM_pn_log_minus_one"=> PhenomHM_TIGER_spinless(-1.),
    #"PhenomHM_pn_log_two_half" => PhenomHM_TIGER_spinless(log(2.5)),
    #"PhenomHM_pn_three"        => PhenomHM_TIGER_spinless(3.),
    "PhenomHM_pn_three_half"   => PhenomHM_TIGER_spinless(3.5)
    )

### calculate fishers ##########################################################
for wf_name in keys(waveforms)
    
    println("\nProcessing: "*wf_name)
    wf = waveforms[wf_name]
    
    folder_name = joinpath(output_folder_name, "resutlts_b_fmin=$(fmin)_fmax=$(fmax)")

    println("\nStoring results into:")
    println(folder_name)
    mkpath(folder_name)

    # extract parameter
    mc     = parameters_mat[1,1:n_events]
    η      = parameters_mat[2,1:n_events]
    χ_1    = parameters_mat[3,1:n_events]
    χ_2    = parameters_mat[4,1:n_events]
    dL     = parameters_mat[5,1:n_events]
    θ      = parameters_mat[6,1:n_events]
    ϕ      = parameters_mat[7,1:n_events]
    iota   = parameters_mat[8,1:n_events]
    ψ      = parameters_mat[9,1:n_events]
    tcoal  = parameters_mat[10,1:n_events]
    Φ_coal = parameters_mat[11,1:n_events]
    z      = parameters_mat[14,1:n_events] 

    # redshift_cutoff
    idx_refshift = z .< 1.0


    println("Calculating Fisher matrices and SNRs")
    @time fisher_matrices, snrs = FisherMatrix(
        wf,
        network,
        mc[idx_refshift], 
        η[idx_refshift], 
        χ_1[idx_refshift], 
        χ_2[idx_refshift], 
        dL[idx_refshift], 
        θ[idx_refshift], 
        ϕ[idx_refshift], 
        iota[idx_refshift], 
        ψ[idx_refshift], 
        tcoal[idx_refshift], 
        Φ_coal[idx_refshift], 
        zeros(Float64, size(mc)), 
        auto_save=false, 
        return_SNR=true, 
        useEarthMotion=true,
        fmin = Float64(fmin),
        fmax = Float64(fmax)
    )
        
    ### postprocessing #########################################################

    # invert fisher matrices and calculate errors 
    covariance_matrices = CovMatrix(fisher_matrices)

    deltak_pn = Array{Float64}(undef, n_events)
    fisher_inverted = Array{Bool}(undef, n_events)
    not_inverted = 0

    for idx_event = 1:n_events
        cmat = covariance_matrices[idx_event, :,:]
        if all(cmat .== 0.0)

            deltak_pn[idx_event] = NaN
            fisher_inverted[idx_event] = false
            not_inverted += 1
        else
            # take the squareroot of the 12th parameter. This is the GR deviation.
            deltak_pn[idx_event] = sqrt(cmat[12, 12])
            fisher_inverted[idx_event] = true
        end
    end

    println("$(not_inverted) Fisher matrices (out of $(n_events)) could not be inverted, or had a SNR below the threshold.")

    ### saving the data #########################################################

    # save everything to .h5 file 
    filename = joinpath(folder_name, "$(wf_name).h5")
    println("Storing: $(filename)")
    h5open(filename, "w") do file
        grp = create_group(file, "parameter")
        write(grp, "mc"      , mc[idx_refshift])
        write(grp, "eta"     , η[idx_refshift])
        write(grp, "chi_1"   , χ_1[idx_refshift])
        write(grp, "chi_2"   , χ_2[idx_refshift])
        write(grp, "dL"      , dL[idx_refshift])
        write(grp, "theta"   , θ[idx_refshift])
        write(grp, "phi"     , ϕ[idx_refshift])
        write(grp, "psi"     , ψ[idx_refshift])
        write(grp, "iota"    , iota[idx_refshift])
        write(grp, "tcoal"   , tcoal[idx_refshift])
        write(grp, "phi_coal", Φ_coal[idx_refshift])
        write(grp, "z"       , z[idx_refshift])

        grp = create_group(file, "results")
        write(grp, "delta_k", deltak_pn)
        write(grp, "snr"    , snrs)
    end

end
