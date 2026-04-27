using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using BenchmarkTools
using MeanFieldGraph
using KrylovKit

## parameters and model
r₊ = 0.5
β = 0.5
λ = 0.5
p = 0.5
model = MarkovChainModel(β * λ, λ, p)

## Dimensions, time duration and data
N = 100
T = 10000
excitatory = MeanFieldGraph.N2excitatory(N, r₊)
data = rand(model, excitatory, T)

## functions to benchmark
aggregated_estimation(data::DiscreteTimeData) = MeanFieldGraph.covariance_vector(data)
function spectral_estimation(data::DiscreteTimeData)
    Σ̂ = MeanFieldGraph.covariance_matrix(data)
    return svdsolve(Σ̂)
end

## benchmark
@benchmark aggregated_estimation($data)
@benchmark spectral_estimation($data)