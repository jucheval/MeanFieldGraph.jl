@testset "model.jl" begin
    @test size(MarkovChainConnectivity(MarkovChainModel(0,0,0), ones(Int, (4,4)))) == 4
    
    data = DiscreteTimeData([true false true; false true true])
    @test length(data) == 3
    @test size(data) == (2,3)
    @test data[1:2].X == [true false; false true]

    @test mvw(MarkovChainModel(1/8,1/2,1/2),1/2) == (1/4, 5/256, 51/256)
end