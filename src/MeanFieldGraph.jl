module MeanFieldGraph

using Clustering: Clustering, KmeansResult, assignments, counts, kmeans
using Distributions: Distributions, Bernoulli, DiscreteUniform, fit, mean
using LinearAlgebra: LinearAlgebra, I, transpose
using Plots: Plots, heatmap, palette, plot

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
