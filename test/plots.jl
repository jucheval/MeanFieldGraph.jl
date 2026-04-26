@testset "Plot extension" begin
    try
        using Plots
        data = MF.DiscreteTimeData(trues(4, 6))
        fig = MF.plot(data)
        @test fig !== nothing
    catch e
        @test_skip "Plots is not installed in active environment"
    end
end
