module HierDist

using Distributions
using Roots
using Trapz

export estimate0Symmetric90CiGaussian, muStdEff4ProdNormal
export hyperparamDistTIGER, getDistributionOnGrid, getNaiveDistributionOnGrid

include("hierarchical_distribution.jl")
include("prod_likelihood_distribution.jl")

end
