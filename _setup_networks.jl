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

function getNetwork(network_name)

    available_networks = ["ETS", "network_0_15km", "network_45_15km"] 
    if !(network_name in available_networks)
        throw(ValueError(network_name, "network not available. Use one of these: $(available_networks)"))
    end

    if network_name == "ETS"
        return [ET]
    elseif network_name == "network_0_15km"
        return [ETLS_15km, ETLMR_0_15km]
    elseif netwoek_name == "network_45_15km"
        return [ETLS_15km, ETLMR_45_15km]
    end 
end 