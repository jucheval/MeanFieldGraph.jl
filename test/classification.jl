@testset "classification.jl" begin
    data = DiscreteTimeData([
        [1 0 0]
        [0 1 1]
        [1 0 0]
        [1 1 0]
    ])

    @test begin
        handcomputation = [
            1 - 2 * 1 / 3, 1 / 2 - 2 * 2 / 3, 1 - 2 * 1 / 3, 3 / 2 - 2 * 2 / 3
        ]
        MF.covariance_vector(data) == handcomputation
    end

    @test begin
        handcomputation_Xs =
            1 / 2 * [
                [0 + 0 0 + 0 0 + 0 0 + 0]
                [1 + 0 0 + 1 1 + 0 1 + 1]
                [0 + 0 0 + 0 0 + 0 0 + 0]
                [1 + 0 0 + 0 1 + 0 1 + 0]
            ]
        handcomputation_Zs =
            1 / 9 * [
                [1 * 1 1 * 2 1 * 1 1 * 2]
                [2 * 1 2 * 2 2 * 1 2 * 2]
                [1 * 1 1 * 2 1 * 1 1 * 2]
                [2 * 1 2 * 2 2 * 1 2 * 2]
            ]
        MF.covariance_matrix(data) == handcomputation_Xs .- handcomputation_Zs
    end

    @test begin
        X = reshape([0.1, 0.2, 0.3, 2.1, 2.2, 2.3], (1, 6))
        R = kmeans(X, 2; init=[1, 6])
        MF.cluster2bool(R) == [false, false, false, true, true, true]
    end

    @test begin
        classification(data) == [true, false, true, true]
    end

    @testset "covariance matrix and vector" begin
        for _ in 1:10
            Tsimu = Int(1e3)
            N = 100
            r₊ = rand()
            excitatory = MF.N2excitatory(N, r₊)
            β = rand()
            λ = 0.3 + 0.7 * rand()
            p = rand()
            model = MarkovChainModel(β * λ, λ, p)
            data = rand(model, N, r₊, Tsimu)
            cov_vec = MeanFieldGraph.covariance_vector(data)
            cov_mat = MeanFieldGraph.covariance_matrix(data)
            @test cov_vec ≈ sum(cov_mat; dims=1)[1, :]
        end
    end
end
