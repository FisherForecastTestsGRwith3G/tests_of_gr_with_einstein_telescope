# Script that defines functions for GR deviations. 
# Please use this script to document the different ways we generate different 
# GR-deviations on top of a catalog. 
#
# There is no import statements in this script
# Inlcude this script into the other script that runs the analysis/simulation 
# for a catalog i.e. use
#    include("_setup_gr_deviations.jl")
# in any other script. 
#

### Comment on functions that create GR parameter
# All functions in this script should take a vectors of the standard gr_parameter,
# as they are used by GW.jl, as inputs. 
#
# The output of these scripts should be a tuple of GR-deviations in the 
# post newtownian coefficients of different order
# delta_pn = (pn-1, pn0, pn0.5, pn1, pn1.5, pn2, log(2.5), 3, log(3.), 3.5)
#
###

"""
Generates the GR-deviations by drawing them with a nomral distribution. 
"""
function deltaPnNormal(
    mc::Vector{Float64}, 
    η::Vector{Float64},
    χ_1::Vector{Float64},
    χ_2::Vector{Float64},
    dL::Vector{Float64},
    θ::Vector{Float64},
    ϕ::Vector{Float64},
    iota::Vector{Float64},
    ψ::Vector{Float64},
    tcoal::Vector{Float64},
    Φ_coal::Vector{Float64},
    lambda1::Union{Vector{Float64},Nothing}=nothing,
    lambda2::Union{Vector{Float64},Nothing}=nothing;
    mu::Union{Vector{Float64},Nothing}=nothing,
    sigma::Union{Vector{Float64},Nothing}=nothing
    )

    mu = mu*ones(Float64, 10)
    sigma = sigma*ones(Float64, 10)
    return deltaPnNormal(mc, η, χ_1, χ_2, dL, θ, ϕ, iota, ψ, tcoal, Φ_coal, lambda1, lambda2, mu=mu, sigma=sigma)
end

function deltaPnNormal(
    mc::Vector{Float64}, 
    η::Vector{Float64},
    χ_1::Vector{Float64},
    χ_2::Vector{Float64},
    dL::Vector{Float64},
    θ::Vector{Float64},
    ϕ::Vector{Float64},
    iota::Vector{Float64},
    ψ::Vector{Float64},
    tcoal::Vector{Float64},
    Φ_coal::Vector{Float64},
    lambda1::Union{Vector{Float64},Nothing}=nothing,
    lambda2::Union{Vector{Float64},Nothing}=nothing;
    mu::Union{Vector{Float64},Nothing}=nothing,
    sigma::Union{Vector{Float64},Nothing}=nothing
    )

    n_events = length(mc)
    # TODO: add check that all inputs have the same length.

    if mu === nothing
        mu = zeros(10)
    end
    if sigma === nothing
        sigma = zeros(10)
    end
    
    delta_pn = [
        mu[1] .+ sigma[1].*randn(n_events),   # pn = -1
        mu[2] .+ sigma[2].*randn(n_events),   # pn = 0
        mu[3] .+ sigma[3].*randn(n_events),   # pn = 0.5
        mu[4] .+ sigma[4].*randn(n_events),    # pn = 1
        mu[5] .+ sigma[5].*randn(n_events),   # pn = 1.5
        mu[6] .+ sigma[6].*randn(n_events),   # pn = 2
        mu[7] .+ sigma[7].*randn(n_events),   # pn = log(2.5)
        mu[8] .+ sigma[8].*randn(n_events),   # pn = 3
        mu[9] .+ sigma[9].*randn(n_events),   # pn = log(3.)
        mu[10] .+ sigma[10].*randn(n_events)  # pn = 3.5
    ]
    
    return delta_pn
end 

"""
Generates the GR-deviations as in farr et al. 
see: 
"""
function farrEtAl(
    mc::Vector{Float64}, 
    η::Vector{Float64},
    χ_1::Vector{Float64},
    χ_2::Vector{Float64},
    dL::Vector{Float64},
    θ::Vector{Float64},
    ϕ::Vector{Float64},
    iota::Vector{Float64},
    ψ::Vector{Float64},
    lambda1::Union{Vector{Float64},Nothing}=nothing,
    lambda2::Union{Vector{Float64},Nothing}=nothing;
    mu::Union{Vector{Float64},Nothing}=nothing,
    sigma::Union{Vector{Float64},Nothing}=nothing 
    )

    n_events = length(η)
    delta_pn0 = 1.1 .- 2.0 .*η
    delta_pn1 = 0.1*ones(n_events)

    delta_pn = (
        zeros(n_events),  # pn = -1
        delta_pn0,        # pn = 0
        zeros(n_events),  # pn = 0.5
        delta_pn1,        # pn = 1
        zeros(n_events),  # pn = 1.5
        zeros(n_events),  # pn = 2
        zeros(n_events),  # pn = log(2.5)
        zeros(n_events),  # pn = 3
        zeros(n_events),  # pn = log(3.)
        zeros(n_events)   # pn = 3.5
    )

    return delta_pn
     
end 