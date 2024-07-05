module MeanFieldGraph

using Distributions, LinearAlgebra, Plots

export MarkovChainModel, DiscreteTimeData, MarkovChainConnectivity
export mvw, mvw_inf
export rand
export estimators, fit
export plot


include("model.jl")
include("simulate.jl")
include("estimation.jl")
include("plots.jl")

end
