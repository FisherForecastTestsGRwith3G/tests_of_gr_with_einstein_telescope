# Script that sets up all the networks for investigations with GW.jl
# Please use this script to document the different detector/network setups we use
#
# There is no import statements in this script
# Inlcude this script into the other script that runs the analysis/simulation 
# for a catalog i.e. use
#    include("_setup_netwoks.jl")
# in any other script. 
#
# TODO: Add better documentation for the networks.
# TODO: Clean up namespace

detector_data_dir = string(parentdir)*"/GW.jl/useful_files/ET_curves/"
detector_data_dir_LVK = string(parentdir)*"/GW.jl/useful_files/"

f_ET_10_1, PSD_ET_10_1 = _readPSD(detector_data_dir*"ET10km.txt", cols=[1, 2])
f_ET_10_2, PSD_ET_10_2 = _readPSD(detector_data_dir*"ET10km.txt", cols=[1, 3])
f_ET_10, PSD_ET_10 = _readPSD(detector_data_dir*"ET10km.txt", cols=[1, 4])

f_ET_15_1, PSD_ET_15_1 = _readPSD(detector_data_dir*"ET15km.txt", cols=[1, 2])
f_ET_15_2, PSD_ET_15_2 = _readPSD(detector_data_dir*"ET15km.txt", cols=[1, 3])
f_ET_15, PSD_ET_15 = _readPSD(detector_data_dir*"ET15km.txt", cols=[1, 4])

f_ET_20_1, PSD_ET_20_1 = _readPSD(detector_data_dir*"ET20km.txt", cols=[1, 2])
f_ET_20_2, PSD_ET_20_2 = _readPSD(detector_data_dir*"ET20km.txt", cols=[1, 3])
f_ET_20, PSD_ET_20 = _readPSD(detector_data_dir*"ET20km.txt", cols=[1, 4])

ET = deepcopy(ETS)
ET.fNoise = f_ET_10
ET.psd = PSD_ET_10

ETLS_10km = deepcopy(ETLS)
ETLS_10km.label = "ETLS_10km"
ETLS_10km.fNoise = f_ET_10
ETLS_10km.psd = PSD_ET_10

ETLMR_0_10km = deepcopy(ETLMR)
ETLMR_0_10km.label = "ETLMR_0_10km"
ETLMR_0_10km.fNoise = f_ET_10
ETLMR_0_10km.psd = PSD_ET_10

orientation_CBC = _orientationBigCircle(ETLS_10km, ETLMR_0_10km)
ETLMR_0_10km.orientation_rad = orientation_CBC - pi/4

ETLMR_45_10km = deepcopy(ETLMR)
ETLMR_45_10km.fNoise = f_ET_10
ETLMR_45_10km.psd = PSD_ET_10
ETLMR_45_10km.orientation_rad = orientation_CBC
ETLMR_45_10km.label = "ETLMR_45_10km"

ETLS_15km = deepcopy(ETLS)
ETLS_15km.label = "ETLS_15km"
ETLS_15km.fNoise = f_ET_15
ETLS_15km.psd = PSD_ET_15

ETLMR_0_15km = deepcopy(ETLMR)
ETLMR_0_15km.label = "ETLMR_0_15km"
ETLMR_0_15km.fNoise = f_ET_15
ETLMR_0_15km.psd = PSD_ET_15
ETLMR_0_15km.orientation_rad = orientation_CBC - pi/4

ETLMR_45_15km = deepcopy(ETLMR)
ETLMR_45_15km.label = "ETLMR_45_15km"
ETLMR_45_15km.fNoise = f_ET_15
ETLMR_45_15km.psd = PSD_ET_15
ETLMR_45_15km.orientation_rad = orientation_CBC

#Setting up LVK (LHV) network, using O3 sensitivities
#Loading realistic O3 PSD for LVK
LIGO_L_O3 = Detector(getCoords(LIGO_L)..., 'L', _readASD(detector_data_dir_LVK * "LVC_O1O2O3/O3-L1-C01_CLEAN_SUB60HZ-1240573680.0_sensitivity_strain_asd.txt")...,  "LIGO_L_O3")
LIGO_H_O3 = Detector(getCoords(LIGO_H)..., 'L', _readASD(detector_data_dir_LVK * "LVC_O1O2O3/O3-H1-C01_CLEAN_SUB60HZ-1251752040.0_sensitivity_strain_asd.txt")...,  "LIGO_H_O3")
VIRGO_O3 = Detector(getCoords(VIRGO)..., 'L', _readASD(detector_data_dir_LVK * "LVC_O1O2O3/O3-V1_sensitivity_strain_asd.txt")...,  "VIRGO_O3")


function getNetwork(network_name)

    available_networks = ["ETS", "network_0_15km", "network_45_15km", "LHV", "LHVK", "LHV_O3"]
    if !(network_name in available_networks)
        throw(ArgumentError("$(network_name) network not available. Use one of these: $(available_networks)"))
    end

    if network_name == "ETS"
        return [ET]
    elseif network_name == "network_0_15km"
        return [ETLS_15km, ETLMR_0_15km]
    elseif network_name == "network_45_15km"
        return [ETLS_15km, ETLMR_45_15km]
    elseif network_name == "LHV"
        return [_available_detectors("LIGO_L"), _available_detectors("LIGO_H"), _available_detectors("VIRGO")]
    elseif network_name == "LHVK"
        return [_available_detectors("LIGO_L"), _available_detectors("LIGO_H"), _available_detectors("VIRGO"), _available_detectors("KAGRA")]
    elseif network_name == "LHV_O3"
        return [LIGO_L_O3, LIGO_H_O3, VIRGO_O3]
    end 
end 