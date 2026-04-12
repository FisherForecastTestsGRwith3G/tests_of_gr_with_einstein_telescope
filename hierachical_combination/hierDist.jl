module HierDist

using Distributions
using Roots
using Trapz
using CairoMakie
import Contour as con

export estimate0Symmetric90CiGaussian, muStdEff4ProdNormal
export hyperparamDistTIGER, getDistributionOnGrid, getNaiveDistributionOnGrid
export quantile1dOGD, getEnclosedIsoVolProbOGD
export getContourLevelOGD, getCredibleContourOGD

include("hierarchical_distribution.jl")
include("prod_likelihood_distribution.jl")
include("stat_utils.jl")
include("plotting_utils.jl")

end
