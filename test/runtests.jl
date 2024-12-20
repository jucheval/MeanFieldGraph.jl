using MeanFieldGraph
using Test
using Distributions, LinearAlgebra
using Random
using Clustering
import MeanFieldGraph as MF

@testset "MeanFieldGraph.jl" begin
    include("model.jl")
    include("simulate.jl")
    include("estimation.jl")
    include("classification.jl")
end
