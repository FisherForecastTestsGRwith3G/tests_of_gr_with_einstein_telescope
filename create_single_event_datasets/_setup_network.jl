#-----------------------------------------------------------------------------# 
# ET detector networks                                                        #
#-----------------------------------------------------------------------------#

"""
Sets up the ET in triangular configuration.  
"""
function setupET10kmT(psd_dir)
    
    f_ET_10, PSD_ET_10 = _readPSD(joinpath(psd_dir,"ET10km.txt"), cols=[1, 4])
    
    ET = deepcopy(ETS)
    ET.fNoise = f_ET_10
    ET.psd = PSD_ET_10

    return [ET]
end

"""
Set up ET network with two L-shaped aligned detectors with 15km arms.
"""
function setupET15km0(psd_dir)

    f_ET_15  , PSD_ET_15 = _readPSD(joinpath(psd_dir,"ET15km.txt"), cols=[1, 4])

    orientation_CBC = _orientationBigCircle(ETLS, ETLMR)

    ETLS_15km = deepcopy(ETLS)
    ETLS_15km.label = "ETLS_15km"
    ETLS_15km.fNoise = f_ET_15
    ETLS_15km.psd = PSD_ET_15

    ETLMR_0_15km = deepcopy(ETLMR)
    ETLMR_0_15km.label = "ETLMR_0_15km"
    ETLMR_0_15km.fNoise = f_ET_15
    ETLMR_0_15km.psd = PSD_ET_15
    ETLMR_0_15km.orientation_rad = orientation_CBC - pi/4

    return [ETLS_15km, ETLMR_0_15km]
end

"""
Set up ET network with two L-shaped 45 degrees misaligned detectors with 15km arms.
"""
function setupET15km45(psd_dir)

    f_ET_15  , PSD_ET_15 = _readPSD(joinpath(psd_dir,"ET15km.txt"), cols=[1, 4])

    orientation_CBC = _orientationBigCircle(ETLS, ETLMR)

    ETLS_15km = deepcopy(ETLS)
    ETLS_15km.label = "ETLS_15km"
    ETLS_15km.fNoise = f_ET_15
    ETLS_15km.psd = PSD_ET_15

    ETLMR_45_15km = deepcopy(ETLMR)
    ETLMR_45_15km.label = "ETLMR_0_15km"
    ETLMR_45_15km.fNoise = f_ET_15
    ETLMR_45_15km.psd = PSD_ET_15
    ETLMR_45_15km.orientation_rad = orientation_CBC

    return [ETLS_15km, ETLMR_45_15km]
end

#-----------------------------------------------------------------------------# 
# LIGO-Virgo detector networks                                                #
#-----------------------------------------------------------------------------#

"""
Set up the O3 LIGO-Virgo-KAGRA network. (No KAGRA detector online at O3)
"""
function setupLHV(asd_dir)

    l_asd_path = joinpath(asd_dir, "O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_strain_asd.txt")
    LIGO_L_O3  = Detector(getCoords(LIGO_L)..., 'L', _readASD(l_asd_path)...,  "LIGO_L_O3")
    
    h_asd_path = joinpath(asd_dir, "O3-H1-C01_CLEAN_SUB60HZ-1251752040.0_sensitivity_strain_asd.txt")
    LIGO_H_O3  = Detector(getCoords(LIGO_H)..., 'L', _readASD(h_asd_path)...,  "LIGO_H_O3")

    v_asd_path = joinpath(asd_dir, "O3-V1_sensitivity_strain_asd.txt") 
    VIRGO_O3   = Detector(getCoords(VIRGO)..., 'L', _readASD(v_asd_path)...,  "VIRGO_O3")

    return [LIGO_L_O3, LIGO_H_O3, VIRGO_O3]
end 

#-----------------------------------------------------------------------------# 
# Select detector network                                                     #
#-----------------------------------------------------------------------------#

"""
Function that selects one of the available detector networks
"""
function getNetwork(network_name)

    psd_data_dir_et  = joinpath(@__DIR__, "psd_data", "et_curves")
    asd_data_dir_lvk = joinpath(@__DIR__, "psd_data", "hlv_curves")

    available_networks = ["ETS", "network_0_15km", "network_45_15km", "LHV_O3"]

    if network_name == "ETS"
        return setupET10kmT(psd_data_dir_et)
    elseif network_name == "network_0_15km"
        return setupET15km0(psd_data_dir_et)
    elseif network_name == "network_45_15km"
        return setupET15km45(psd_data_dir_et)
    elseif network_name == "LHV"
        return setupLHV(asd_data_dir_lvk)
    else 
        throw(ArgumentError("$(network_name) network not available. Use one of these: $(available_networks)"))
    end 
end 
