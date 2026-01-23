using Aqua
using Clustering
using Distributions
using Documenter
using JET
using LinearAlgebra
using MeanFieldGraph
import MeanFieldGraph as MF
using Random
using Test

@testset verbose = true "MeanFieldGraph.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(MeanFieldGraph; ambiguities=false, deps_compat=(; check_extras=false))
    end
    @testset "Type stability (JET.jl)" begin
        # test on 1.11 at a minimum and not on pre-release 
        if VERSION >= v"1.12" && isempty(VERSION.prerelease)
            JET.test_package(MeanFieldGraph; target_modules=(MeanFieldGraph,))
        end
    end
    @testset "Doctests" begin
        doctest(MeanFieldGraph)
    end
    include("model.jl")
    include("simulate.jl")
    include("estimation.jl")
    include("classification.jl")
end
