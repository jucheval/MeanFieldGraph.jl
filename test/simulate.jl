@testset "simulate.jl" begin
    @test begin
        testvalue = [true, false, true, false]
        MF.forward_simulation!(
            testvalue,
            MF.MarkovChainConnectivity(MarkovChainModel(0.5, 0.5, 0.0), ones(Int, (4, 4))),
            [true, false, true, false],
        )
        testvalue == [true, true, true, true]
    end

    @test begin
        testvalue = [false, true, false, true]
        MF.forward_simulation!(
            testvalue,
            MF.MarkovChainConnectivity(MarkovChainModel(0, 0.5, 0.0), ones(Int, (4, 4))),
            [true, false, true, false],
        )
        testvalue == [false, false, false, false]
    end

    @test begin
        Random.seed!(1)
        Nsimu = 1e6
        β = rand()
        λ = 0.3 + 0.7 * rand()
        μ = β * λ
        θ = rand(MF.ErdosRenyiGraph(4, rand()))
        excitatory = [true, true, false, false]
        model = MF.MarkovChainConnectivity(MarkovChainModel(μ, λ, 0.0), θ)

        mhat = zeros(4)
        for i in 1:Nsimu
            mhat += MF.stationary_initial_condition(model, excitatory) / Nsimu
        end

        A = 1 / 4 * (θ .* transpose(-1 .+ 2 * excitatory))
        Q = inv(I - (1 - λ) * A)
        L⁻ = sum(A[:, 3:4]; dims=2)
        mtheo = μ * Q * ones(4) - (1 - λ) * Q * L⁻

        sum(abs.(mhat - mtheo)) < 5e-2
    end

    @test begin
        Random.seed!(1)
        Tsimu = Int(1e6)
        β = rand()
        λ = 0.2 + 0.8 * rand()
        μ = β * λ
        θ = rand(MF.ErdosRenyiGraph(4, rand()))
        excitatory = [true, true, false, false]
        model = MF.MarkovChainConnectivity(MarkovChainModel(μ, λ, 0.0), θ)

        data = rand(model, excitatory, Tsimu)
        mhat = mean(data.X; dims=2)

        A = 1 / 4 * (θ .* transpose(-1 .+ 2 * excitatory))
        Q = inv(I - (1 - λ) * A)
        L⁻ = sum(A[:, 3:4]; dims=2)
        mtheo = μ * Q * ones(4) - (1 - λ) * Q * L⁻

        sum(abs.(mhat - mtheo)) < 3e-3
    end
end
