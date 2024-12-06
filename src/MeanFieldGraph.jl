module MeanFieldGraph

using Distributions, LinearAlgebra, Plots, Clustering

export MarkovChainModel, DiscreteTimeData, MarkovChainConnectivity
export mvw, mvw_inf
export rand
export estimators, fit
export classification
export plot

include("model.jl")
include("simulate.jl")
include("estimation.jl")
include("classification.jl")
include("plots.jl")

end
