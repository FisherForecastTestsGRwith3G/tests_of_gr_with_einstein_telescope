"""
DESCRIPTION: 

This module provides functionalities to create the fisher approximation
of the posteriors distributions for every event in a simulated catalog 
of Binary Black Holes BBH. 
"""
module createSED

# -------------------------------------------------------------------------- #
# Preamble                                                                   #
# -------------------------------------------------------------------------- #
using GW 
using Random
using LaTeXStrings

import Base: length

# define all auxiliary functions for the script 
include("./_setup_network.jl")
include("./_setup_catalog.jl")

export BBHCatalog
export applyRedshiftCut, truncateCatalog, merge, createBGRDeviations
export createSEDfromCatalog
export pnoString, pnoNum, pnoLatex, pno_list 

const pno_dict = Dict(
    "-1"       => (-1.0    ,"minus_one"    , L"\delta \phi_{-2}"  )     ,
    "0"        => (0.0     ,"zero"         , L"\delta \phi_{0}"   )     , 
    "0.5"      => (0.5     ,"half"         , L"\delta \phi_{1}" )     ,
    "1"        => (1.0     ,"one"          , L"\delta \phi_{2}"   )     ,   
    "1.5"      => (1.5     ,"one_half"     , L"\delta \phi_{3}" )     , 
    "2"        => (2.0     ,"two"          , L"\delta \phi_{4}"   )     ,
    "log(2.5)" => (log(2.5),"log_two_half" , L"\delta \phi_{5\ell}" ),
    "3"        => (3.0     ,"three"        , L"\delta \phi_{6}"   )     ,
    "log(3.)"  => (log(3.) ,"log_three"    , L"\delta \phi_{6\ell}")  ,
    "3.5"      => (3.5     ,"three_half"   , L"\delta \phi_{7}" )     ,
)
const pno_list = [ "-1","0","0.5","1","1.5","2","log(2.5)","3","log(3.)","3.5"]

function pnoString(pno)
    return pno_dict[pno][2]
end

function pnoNum(pno)
    return pno_dict[pno][1]
end

function pnoLatex(pno)
    return pno_dict[pno][3]
end

# -------------------------------------------------------------------------- #
# define main functions                                                      #
# -------------------------------------------------------------------------- #

"""
INPUTS
--------
    catalog         Array of GR parameters of Binary Black Holes
    pno             Post newtonian order
    wf_family       Waveform family: PhenomD, PhenomHM
    network         Detector network
    mu              Gaussian deviation population parameter: Mean
    sigma           Gaussian deviation population parameter: Standard deviation

OUTPUTS
--------
    fisher::Vector{Matrix{Float64}}
        List of fisher matrices for each event.                 
    
    snr::Vector{Float64}            
        SNR of the signal 
    
    isnr::Vector{Float64}
        Inspiral SNR of the signal   
    
    invc::Vector{Bool}
        Invertibility check of fisher 
    
    delta_k::Vector{Float64}
        standard deviation of fisher approximated single event deviation posterior

    dphi_k::Vector{Float}
        mean value of fisher approximated single event deviation posterior
"""
function createSEDfromCatalog(
    catalog::BBHCatalog,
    network_name::String,
    pno::String,
    pn_deviation::Vector{Float64}, 
    wf_family::String, 
    fmin::Float64,
    seed::Int64
    )

    n_events = length(catalog)
    network  = getNetwork(network_name)
    
    if wf_family == "PhenomD"
        wf_model = PhenomD_TIGER_spinless(pnoNum(pno))
    elseif wf_family == "PhenomHM"
        wf_model = PhenomHM_TIGER_spinless(pnoNum(pno))
    else
        error("Waveform family is not available: $(wf_family)")
    end

    # --------------------------------------------------#
    # Fisher matrices and SNR                           #
    # --------------------------------------------------#
    println("Calculating Fisher matrices and SNRs")
    @time fisher, snr = FisherMatrix(
        wf_model              ,
        network               ,
        catalog.mc            , 
        catalog.eta           , 
        catalog.chi_1         , 
        catalog.chi_2         ,  
        catalog.dL            , 
        catalog.theta         , 
        catalog.phi           , 
        catalog.iota          , 
        catalog.psi           , 
        catalog.t_coal        , 
        catalog.phi_coal      , 
        pn_deviation          , 
        auto_save      = false, 
        return_SNR     = true , 
        useEarthMotion = true ,
        fmin           = fmin
    )

    # --------------------------------------------------#
    # Inspiral SNR                                      #
    # --------------------------------------------------#
    println("Calculating inspiral SNRs")
    # Evaluate SNR only for the Inspiral part of the waveform
    # Truncate waveform in frequency space. For the definition 
    # of the inspiral phase see the paper for the Phenom waveform 
    # models: 
    f_inspiral_cutoff = @. 0.018 / (  catalog.mc / catalog.eta^(3. /5.) ) / GMsun_over_c3

    # Ensure that f_inspiral_cutoff is above fmin
    f_inspiral_cutoff = max.(f_inspiral_cutoff, fmin)

    @time isnr = SNR(
        wf_model                           ,
        network                            ,
        catalog.mc                         , 
        catalog.eta                        , 
        catalog.chi_1                      , 
        catalog.chi_2                      ,  
        catalog.dL                         , 
        catalog.theta                      , 
        catalog.phi                        , 
        catalog.iota                       , 
        catalog.psi                        , 
        catalog.t_coal                     , 
        pn_deviation                       , 
        auto_save       = false            , 
        useEarthMotion  = true             ,
        fmax            = f_inspiral_cutoff,
        fmin            = fmin
    )    

    # --------------------------------------------------#
    # Inversion of Fisher matrices and estimate delta_k #
    # --------------------------------------------------#
    cov_mat         = CovMatrix(fisher)
    delta_k         = Array{Float64}(undef, n_events)
    invc            = Array{Bool}(undef, n_events)
    n_not_inverted  = 0

    for idx_event = 1:n_events
        cm = cov_mat[idx_event, :,:]
        if all(cm .== 0.0)
            # TODO: Find a better way to deal with non-invertible fisher matrices
            delta_k[idx_event] = NaN
            invc[idx_event]    = false
            n_not_inverted     += 1
        else 
            delta_k[idx_event] = sqrt(cm[12, 12])
            invc[idx_event]    = true
        end
    end

    # --------------------------------------------------#
    # Estimate dphi_k                                   #
    # --------------------------------------------------#
    dphi_k = pn_deviation .+ delta_k .* randn(n_events)

    return fisher, snr, isnr, invc, delta_k, dphi_k
end

end #module
