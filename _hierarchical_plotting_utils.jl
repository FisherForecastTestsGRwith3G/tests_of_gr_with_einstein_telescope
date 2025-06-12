# using Trapz
import Contour as con


"""
Calculate percentiles of d-dimensional discrete probability distribution
by scannining the probabilities from top to bottom and summing the up.  
"""
function calcPercentileLvl(    
    p_values::Array{Float64}, 
    percentiles=[0.3935, 0.8647, 0.9889],
    weights=nothing
    )

    if weights == nothing
        weights = ones(size(p_values))
    end

    if any(percentiles .>= 1) | any(percentiles .<= 0)
         throw(ValueError("Percentiles must be values between 0 and 1!"))
    end

    #calculate contours 
    bin_values = p_values[:] .* weights[:]
    bin_values = sort(bin_values)
    z = cumsum(bin_values)
    z = z / z[end]

    percentile_level = []
    percentiles = sort(percentiles, rev=true)
    for perc in percentiles
        idx = findfirst(z .>= (1.0-perc))

        if !(idx == nothing)
            append!(percentile_level, [bin_values[idx]])
        else
            throw(ErrorException("Contour level could not be found. Have a look at the posterior handed!"))
        end
    end

    return percentile_level
end

"""
Summarizes the reconstructed hierarchical distribution in a single plot.
"""
function distributionSummaryPlot(
    mu_values,
    sig_values,
    p_mu_sig,
    p_sig, 
    p_mu; 
    val_inj=nothing,
    title=""
    ) 

    # calculate limits
    mu_lim = [minimum(mu_values), maximum(mu_values)]
    sig_lim = [minimum(sig_values), maximum(sig_values)]
    pmu_lim = [0, maximum(p_mu)*1.2]
    psig_lim = [0, maximum(p_sig)*1.2]

    # calculate contours 
    percentiles = [0.3935, 0.8647, 0.9889]
    percentile_level = calcPercentileLvl(p_mu_sig, percentiles)

    # draw contour plot
    percentile_level = vcat(percentile_level, [maximum(p_mu_sig)])
    percentiles = vcat(percentiles, [1.])
    cplot = contour(
        mu_values,
        sig_values,
        transpose(p_mu_sig),
        fill=true,
        color=:Blues_3,
        cbar = false, 
        levels = percentile_level
    )
    if !isnothing(val_inj)
        plot!(cplot, val_inj[1]*[1,1], sig_lim, color=:orange, linewidth = 2, label = "injected value")
        plot!(cplot, mu_lim, val_inj[2]*[1,1], color=:orange, linewidth = 2, label = "")
    end
    xlims!(cplot, mu_lim...)
    ylims!(cplot, sig_lim...)
    ylabel!(cplot,"σ")
    title!(cplot, title)
    
    sig_plot = plot(p_sig, sig_values, label = "", ticks = :native, ytickfontcolor = RGBA(0,0,0,0))
    if !isnothing(val_inj)
        plot!(sig_plot, psig_lim, val_inj[2]*[1,1], color=:orange, linewidth = 2,label = "")
    end
    xlims!(sig_plot, psig_lim...)
    ylims!(sig_plot, sig_lim...)
    xlabel!(sig_plot, "p(σ|D)")
    
    mu_plot = plot(mu_values, p_mu, label = "")
    if !isnothing(val_inj)
        plot!(mu_plot, val_inj[1]*[1,1], pmu_lim, color=:orange, linewidth = 2,label = "")
    end
    xlims!(mu_plot, mu_lim...)
    ylims!(mu_plot, pmu_lim...)
    xlabel!(mu_plot, "μ")
    ylabel!(mu_plot, "p(μ|D)")

    l = @layout [ jeff{0.7w, 0.8h} karl{0.3w, 0.8h} ; hubert{0.7w, 0.2h} gianlu{0.3w, 0.2h}]
     
    return plot(cplot, sig_plot, mu_plot, layout = l)

end

"""
Calculates a line in the mu-sigma plane, representing the 90CI contour
"""
function calculate90CIContour(mu_values, sig_values, p_mu_sig) 

    # calculate the probability that corresponds to the 90CI
    c_levels = calcPercentileLvl(p_mu_sig, [0.90])
    
    # determine the 90CI isoline
    cl_x = []
    cl_y = []
    for cl in con.levels(con.contours(mu_values, sig_values, p_mu_sig, c_levels))
        lvl = con.level(cl) # the z-value of this contour level
        for line in con.lines(cl)
            xs, ys = con.coordinates(line) # coordinates of this line segment
            append!(cl_x, [xs])
            append!(cl_y, [ys])
        end
    end

    return cl_x, cl_y
end