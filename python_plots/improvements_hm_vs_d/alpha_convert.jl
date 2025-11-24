using Pkg
using HDF5
import JSON
using QuadGK

################################################################################
### auxilary funcion                                                         ###
################################################################################
function get_H_z(z, H0, Omega0_m, Omega0_Lambda)
    H_z = H0 * (Omega0_m * (1 .+ z) .^ 3 .+ Omega0_Lambda) .^ 0.5 # km/s 1/Mpc
    return H_z
end

function get_dL(z)
        
    clight= 2.99792458e5 # km/s
    H0 = 67.66 # km/s 1/Mpc
    Omega0_m = 0.3153
    Omega0_Lambda = 1 - Omega0_m
    dL = zeros(length(z))
    jj = 1
    for ii in z
        dL[jj] =
            quadgk(k -> 1.0 ./ get_H_z(k, H0, Omega0_m, Omega0_Lambda), 0.0, ii)[1] *
            clight *
            (1.0 + ii)
        jj += 1
    end
    return dL
end

function get_z(dL)
    # use bisection method to find z
    z=1
    zmin = 0.
    zmax = 30.
    dL = dL*1e3 # convert to Mpc
    iteration = 0
    while abs.(dL - get_dL(z)[1]) .> 1e-3
        iteration += 1
        if iteration > 100
            println("Could not find z")
            break
        end
        z = 0.5*(zmin + zmax)
        if dL .> get_dL(z)[1]
            zmin = z
        else
            zmax = z
        end
    end
    return z
end 

################################################################################
### Load all the data                                                        ###
################################################################################

data_folder_name = Dict()
data_folder_name["D"] = "/mnt/c/Users/Utente/Desktop/GitRepositories/BGR_with_GWJulia/supplementary/fmin_investigation/paper_samples/final_run_LVK_PhenomD/data/LHV_O3"
data_folder_name["HM"]  = "/mnt/c/Users/Utente/Desktop/GitRepositories/BGR_with_GWJulia/supplementary/fmin_investigation/paper_samples/final_run_LVK/data/LHV_O3"

### get global indices 
global_index = Dict()
data = Dict()
snr_data = Dict()
isnr_data = Dict()

pno_list = ["pn_zero", "pn_one", "pn_three", "pn_three_half", "pn_log_two_half"]
waveforms = ["D", "HM"]

for pno in pno_list
    for wf in waveforms

        file_name = joinpath(data_folder_name[wf] , pno, "single_event_measurements.h5")
        data[wf] = Dict()
        h5open(file_name, "r") do data_file
            data[wf]["delta_k"] = read(data_file, "delta_k")
        end

        file_name = joinpath(data_folder_name[wf] , pno, "snrs.h5")
        snr_data[wf] = Dict()
        h5open(file_name, "r") do snr_file
            snr_data[wf]["values"] = read(snr_file, "values")  
        end

        file_name = joinpath(data_folder_name[wf] , pno, "inspiral_snrs.h5")
        isnr_data[wf] = Dict()
        h5open(file_name, "r") do snr_file
            isnr_data[wf]["values"] = read(snr_file, "values")  
        end

        catalog = Dict()
        file_name = "/mnt/c/Users/Utente/Desktop/GitRepositories/BGR_with_GWJulia/supplementary/fmin_investigation/paper_samples/final_run_LVK/catalog_w_deviations.h5"
        h5open(file_name, "r") do catalog_file
            grp = catalog_file["parameter"]
            catalog["mc"]            = read(grp, "mc")  
            catalog["eta"]           = read(grp, "eta")  
            catalog["chi1"]          = read(grp, "chi1")  
            catalog["chi2"]          = read(grp, "chi2")  
            catalog["dL"]            = read(grp, "dL")  
            catalog["theta"]         = read(grp, "theta")  
            catalog["phi"]           = read(grp, "phi")  
            catalog["iota"]          = read(grp, "iota")  
            catalog["psi"]           = read(grp, "psi")  
        end

        ### Save data

        # save everything to .h5 file 
        wf_name = ""
        if wf == "D"
            wf_name = "PhenomD"
        elseif wf == "HM"
            wf_name = "PhenomHM"
        end
        folder_name = "/mnt/c/Users/Utente/Desktop/GitRepositories/BGR_with_GWJulia/supplementary/fmin_investigation/data/results_paper"
        mkpath(folder_name)
        filename = joinpath(folder_name, "$(wf_name)_$(pno).h5")
        println("Storing: $(filename)")
        h5open(filename, "w") do file
            grp = create_group(file, "parameter")
            write(grp, "mc"   , catalog["mc"])
            write(grp, "eta"  , catalog["eta"])
            write(grp, "chi_1", catalog["chi1"])
            write(grp, "chi_2", catalog["chi2"])
            write(grp, "dL"   , catalog["dL"])
            write(grp, "theta", catalog["theta"])
            write(grp, "phi"  , catalog["phi"])
            write(grp, "psi"  , catalog["psi"] )
            write(grp, "iota" , catalog["iota"])
            write(grp, "z"    , get_z.(catalog["dL"]))

            grp = create_group(file, "results")
            write(grp, "delta_k", data[wf]["delta_k"])
            write(grp, "snr"    , snr_data[wf]["values"])
            write(grp, "isnr"   , isnr_data[wf]["values"])
        end
    end
end

