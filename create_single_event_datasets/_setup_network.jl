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
    ETLMR_45_15km.label = "ETLMR_45_15km"
    ETLMR_45_15km.fNoise = f_ET_15
    ETLMR_45_15km.psd = PSD_ET_15
    ETLMR_45_15km.orientation_rad = orientation_CBC

    return [ETLS_15km, ETLMR_45_15km]
end

#-----------------------------------------------------------------------------# 
# LIGO-Virgo detector networks                                                #
#-----------------------------------------------------------------------------#

function _setupHLVNetwork(l_asd_path, h_asd_path, v_asd_path, label_suffix)

    ligo_l = Detector(getCoords(LIGO_L)..., 'L', _readASD(l_asd_path)..., "LIGO_L_$(label_suffix)")
    ligo_h = Detector(getCoords(LIGO_H)..., 'L', _readASD(h_asd_path)..., "LIGO_H_$(label_suffix)")
    virgo  = Detector(getCoords(VIRGO)...,  'L', _readASD(v_asd_path)..., "VIRGO_$(label_suffix)")

    return [ligo_l, ligo_h, virgo]
end

"""
Set up the O3a HLV network.
"""
function setupHLVO3a(asd_dir)

    l_asd_path = joinpath(asd_dir, "O3a", "O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_strain_asd.txt")
    h_asd_path = joinpath(asd_dir, "O3a", "O3-H1-C01_CLEAN_SUB60HZ-1251752040.0_sensitivity_strain_asd.txt")
    v_asd_path = joinpath(asd_dir, "O3a", "O3-V1_sensitivity_strain_asd.txt")

    return _setupHLVNetwork(l_asd_path, h_asd_path, v_asd_path, "O3a")
end

"""
Set up the O3b HLV network.
"""
function setupHLVO3b(asd_dir)

    l_asd_path = joinpath(asd_dir, "O3b", "O3-L1-C01_CLEAN_SUB60HZ-1262141640.0_sensitivity_strain_asd.txt")
    h_asd_path = joinpath(asd_dir, "O3b", "O3-H1-C01_CLEAN_SUB60HZ-1262197260.0_sensitivity_strain_asd.txt")
    v_asd_path = joinpath(asd_dir, "O3b", "O3-V1-1265246178_sensitivity_strain_asd.txt")

    return _setupHLVNetwork(l_asd_path, h_asd_path, v_asd_path, "O3b")
end

"""
Set up the pre-O4 estimate HLV network used as pO4/O4a sensitivity estimate.
"""
function setupHLVpO4(asd_dir)

    l_asd_path = joinpath(asd_dir, "O4_preO4estimates", "aligo_O4high.txt")
    h_asd_path = joinpath(asd_dir, "O4_preO4estimates", "aligo_O4high.txt")
    v_asd_path = joinpath(asd_dir, "O4_preO4estimates", "avirgo_O4high_NEW.txt")

    return _setupHLVNetwork(l_asd_path, h_asd_path, v_asd_path, "pO4")
end

setupHLVO3(asd_dir) = setupHLVO3a(asd_dir)

#-----------------------------------------------------------------------------# 
# Select detector network                                                     #
#-----------------------------------------------------------------------------#

"""
Function that selects one of the available detector networks
"""
function getNetwork(network_name)

    psd_data_dir_et  = joinpath(@__DIR__, "psd_data", "et_curves")
    asd_data_dir_lvk = joinpath(@__DIR__, "psd_data", "hlv_curves")

    available_networks = ["ETS", "ET_0_15km", "ET_45_15km", "HLV", "HLV_O3", "HLV_O3a", "HLV_O3b", "HLV_pO4"]

    if network_name == "ETS"
        return setupET10kmT(psd_data_dir_et)
    elseif network_name == "ET_0_15km"
        return setupET15km0(psd_data_dir_et)
    elseif network_name == "ET_45_15km"
        return setupET15km45(psd_data_dir_et)
    elseif network_name == "HLV_O3a"
        return setupHLVO3a(asd_data_dir_lvk)
    elseif network_name == "HLV_O3b"
        return setupHLVO3b(asd_data_dir_lvk)
    elseif network_name == "HLV_pO4"
        return setupHLVpO4(asd_data_dir_lvk)
    else 
        throw(ArgumentError("$(network_name) network not available. Use one of these: $(available_networks)"))
    end 
end 
