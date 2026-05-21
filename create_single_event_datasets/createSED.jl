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
export computeSNRsFromCatalog
export pnoString, pnoNum, pnoLatex, pno_list 

const pno_dict = Dict(
    "-1"       => (-1.0    ,"minus_one"    , L"\delta \varphi_{-2}"  )     ,
    "0"        => (0.0     ,"zero"         , L"\delta \varphi_{0}"   )     , 
    "0.5"      => (0.5     ,"half"         , L"\delta \varphi_{1}" )     ,
    "1"        => (1.0     ,"one"          , L"\delta \varphi_{2}"   )     ,   
    "1.5"      => (1.5     ,"one_half"     , L"\delta \varphi_{3}" )     , 
    "2"        => (2.0     ,"two"          , L"\delta \varphi_{4}"   )     ,
    "log(2.5)" => (log(2.5),"log_two_half" , L"\delta \varphi_{5\ell}" ),
    "3"        => (3.0     ,"three"        , L"\delta \varphi_{6}"   )     ,
    "log(3.)"  => (log(3.) ,"log_three"    , L"\delta \varphi_{6\ell}")  ,
    "3.5"      => (3.5     ,"three_half"   , L"\delta \varphi_{7}" )     ,
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

function getWaveformModel(wf_family::String, pno::String)
    if wf_family == "PhenomD"
        return PhenomD_TIGER_spinless(pnoNum(pno))
    elseif wf_family == "PhenomHM"
        return PhenomHM_TIGER_spinless(pnoNum(pno))
    else
        error("Waveform family is not available: $(wf_family)")
    end
end

function computeSNRsFromCatalog(
    catalog::BBHCatalog,
    network_name::String,
    pno::String,
    wf_family::String,
    fmin::Float64,
)
    network = getNetwork(network_name)
    wf_model = getWaveformModel(wf_family, pno)
    zero_pn_deviation = zeros(Float64, length(catalog))

    println("Starting SNR calculations for waveform family $(wf_family) on network $(network_name)")
    println("Calculating network SNRs at zero deviation")
    @time _, snr = FisherMatrix(
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
        zero_pn_deviation     ,
        auto_save      = false,
        return_SNR     = true ,
        useEarthMotion = true ,
        fmin           = fmin
    )

    println("Calculating inspiral-only SNRs")
    f_inspiral_cutoff = @. 0.018 / (catalog.mc / catalog.eta^(3. / 5.)) / GMsun_over_c3
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
        zero_pn_deviation                  ,
        auto_save       = false            ,
        useEarthMotion  = true             ,
        fmax            = f_inspiral_cutoff,
        fmin            = fmin
    )

    return snr, isnr
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
    seed::Int64;
    precomputed_snr::Union{Nothing, Vector{Float64}}=nothing,
    precomputed_isnr::Union{Nothing, Vector{Float64}}=nothing,
    snr_threshold::Float64=0.0,
    inspiral_snr_threshold::Float64=0.0,
    )

    n_events = length(catalog)
    network  = getNetwork(network_name)
    
    wf_model = getWaveformModel(wf_family, pno)
    println("Running single-event Fisher analysis for PN order $(pno) and waveform family $(wf_family)")

    # --------------------------------------------------#
    # Fisher matrices and SNR                           #
    # --------------------------------------------------#
    if isnothing(precomputed_snr)
        println("No precomputed SNRs were provided; calculating Fisher matrices and SNRs for all events")
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
        selected_idx = trues(n_events)

    else
        snr = precomputed_snr
        selected_idx = (snr .> snr_threshold) .& (precomputed_isnr .> inspiral_snr_threshold)

        println("Using precomputed SNRs with thresholds snr > $(snr_threshold) and isnr > $(inspiral_snr_threshold)")
        println("Calculating Fisher matrices for $(sum(selected_idx)) / $(n_events) selected events")
        if any(selected_idx)
            @time fisher_selected = FisherMatrix(
                wf_model                      ,
                network                       ,
                catalog.mc[selected_idx]      , 
                catalog.eta[selected_idx]     , 
                catalog.chi_1[selected_idx]   ,  
                catalog.chi_2[selected_idx]   ,  
                catalog.dL[selected_idx]      , 
                catalog.theta[selected_idx]   , 
                catalog.phi[selected_idx]     , 
                catalog.iota[selected_idx]    , 
                catalog.psi[selected_idx]     ,  
                catalog.t_coal[selected_idx]  , 
                catalog.phi_coal[selected_idx], 
                pn_deviation[selected_idx]    , 
                auto_save      = false        , 
                return_SNR     = false        ,
                useEarthMotion = true         ,
                fmin           = fmin
            )

            n_parameters = size(fisher_selected, 2)
            fisher = zeros(Float64, n_events, n_parameters, n_parameters)
            fisher[selected_idx, :, :] = fisher_selected

        else
            fisher = zeros(Float64, n_events, 12, 12)
        end
    end

    # --------------------------------------------------#
    # Inspiral SNR                                      #
    # --------------------------------------------------#
    if isnothing(precomputed_isnr)
        println("No precomputed inspiral SNRs were provided; calculating them now")
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
    else
        isnr = precomputed_isnr
    end

    # --------------------------------------------------#
    # Inversion of Fisher matrices and estimate delta_k #
    # --------------------------------------------------#
    println("Building covariance matrices for the selected Fisher matrices")
    cov_mat = zeros(Float64, size(fisher))
    cov_mat[selected_idx, :, :] = CovMatrix(fisher[selected_idx, :, :])
    
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
    rng = Random.MersenneTwister(seed)
    dphi_noise = Vector{Float64}(undef, n_events)
    Random.randn!(rng, dphi_noise)
    dphi_k = pn_deviation .+ delta_k .* dphi_noise

    return fisher, snr, isnr, invc, delta_k, dphi_k
end

end #module
