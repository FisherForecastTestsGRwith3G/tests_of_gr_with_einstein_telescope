"""
Catalog of BBH events
"""
struct BBHCatalog
    mc::Array{Float64}       # Chirp mass in detector frame     
    eta::Array{Float64}      # Symmetric mass ratio
    chi_1::Array{Float64}    # Z-component of primary spin
    chi_2::Array{Float64}    # Z-component of secondary spin
    dL::Array{Float64}       # Luminosity distance
    theta::Array{Float64}    #    
    phi::Array{Float64}      # 
    iota::Array{Float64}     # Inclination  
    psi::Array{Float64}      # Polarization angle
    t_coal::Array{Float64}   # Time of coalescence
    phi_coal::Array{Float64} # Coalescence phase
    z::Array{Float64}        # Redshift
end

"""
Removes events with redshift above a threshold.
"""
function applyRedshiftCut(catalog::BBHCatalog, z_thresh::Float64)

    idx_redshift = catalog.z .<= z_thresh

    return BBHCatalog(
        catalog.mc[idx_redshift],
        catalog.eta[idx_redshift],
        catalog.chi_1[idx_redshift],
        catalog.chi_2[idx_redshift],
        catalog.dL[idx_redshift],
        catalog.theta[idx_redshift],
        catalog.phi[idx_redshift],
        catalog.iota[idx_redshift],
        catalog.psi[idx_redshift],
        catalog.t_coal[idx_redshift],
        catalog.phi_coal[idx_redshift],
        catalog.z[idx_redshift]
    )
end

"""
Truncate catalog to only have specified number of events
"""
function truncateCatalog(catalog::BBHCatalog, n_events)
    return BBHCatalog(
        catalog.mc[1:n_events],
        catalog.eta[1:n_events],
        catalog.chi_1[1:n_events],
        catalog.chi_2[1:n_events],
        catalog.dL[1:n_events],
        catalog.theta[1:n_events],
        catalog.phi[1:n_events],
        catalog.iota[1:n_events],
        catalog.psi[1:n_events],
        catalog.t_coal[1:n_events],
        catalog.phi_coal[1:n_events],
        catalog.z[1:n_events]
    )
end

"""
Returns the number of events of the catalog
"""
function length(catalog::BBHCatalog)
    return length(catalog.mc)
end

"""
Exterior constructor for BBHCatalog
"""
function BBHCatalog(n_events, seed)

    gr_catalog = GenerateCatalog(
        n_events,
        "BBH"; 
        time_delay_in_Myr = 10., 
        seed_par = seed, 
        SFR = "Madau&Dickinson", 
        name_catalog = nothing, 
        local_rate = nothing, 
        auto_save=false # not available in the current version of GWInference
        )

    # z_catalog = get_z.(gr_catalog[5])
    return BBHCatalog(
        gr_catalog[1],
        gr_catalog[2],
        gr_catalog[3],
        gr_catalog[4],
        gr_catalog[5],
        gr_catalog[6],
        gr_catalog[7],
        gr_catalog[8],
        gr_catalog[9],
        gr_catalog[10],
        gr_catalog[11],
        gr_catalog[14]
    )
end 

"""
Obtain z of the catalog from dL and cosmology
"""
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

"""
Merge catalogs
"""
function merge(catalog1::BBHCatalog, catalog2::BBHCatalog)
    return BBHCatalog(ntuple(fieldcount(BBHCatalog)) do i
        values1 = getfield(catalog1, i)
        values2 = getfield(catalog2, i)
        merged_values = Vector{Float64}(undef, length(values1) + length(values2))
        copyto!(merged_values, 1, values1, 1, length(values1))
        copyto!(merged_values, length(values1) + 1, values2, 1, length(values2))
        merged_values
    end...)
end

"""
Function that adds a deviation parameter to the catalog.
The deviation follows a Gaussian population distribution, specified by the hyperparameter
mu (mean value of the deviation) and sigma (std of the deviations).  
"""
function createBGRDeviations(catalog::BBHCatalog, mu::Float64, sigma::Float64, seed::Int=1234)
    
    n_events = length(catalog)

    rng = Random.MersenneTwister(seed)
    deviations = Vector{Float64}(undef, n_events)
    Random.randn!(rng, deviations)
    deviations .*= sigma
    deviations .+= mu

    return deviations
end
