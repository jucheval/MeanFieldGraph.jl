@testset "classification.jl" begin
    @test begin
        data = DiscreteTimeData([
            [1 0 0]
            [0 1 1]
            [1 0 0]
            [1 1 0]
        ])
        handcomputation = [
            1 - 2 * 1 / 3, 1 / 2 - 2 * 2 / 3, 1 - 2 * 1 / 3, 3 / 2 - 2 * 2 / 3
        ]
        MF.covariance_vector(data) == handcomputation
    end

    @test begin
        X = reshape([0.1, 0.2, 0.3, 2.1, 2.2, 2.3], (1, 6))
        R = kmeans(X, 2; init=[1, 6])
        MF.cluster2bool(R) == [false, false, false, true, true, true]
    end

    @test begin
        data = DiscreteTimeData([
            [1 0 0]
            [0 1 1]
            [1 0 0]
            [1 1 0]
        ])
        classification(data) == [true, false, true, true]
    end
end
