@testset "estimation.jl" begin
    @test begin
        Random.seed!(1)
        Tsimu = Int(1e4)
        N = 1000
        r₊ = rand()
        excitatory = MF.N2excitatory(N, r₊)
        β = rand()
        λ = .3 + .7*rand()
        p = rand()
        model = MarkovChainModel(β*λ,λ,p)
        theo = mvw(model, r₊)

        hat = (0,0,0)
        hat_inf = (0,0,0)
        Nsimu = 1e2
        for i in 1:Nsimu
            θ = rand(MF.ErdosRenyiGraph(N, p))
            modelconnec = MarkovChainConnectivity(model, θ) 
            data = rand(modelconnec, excitatory, Tsimu)
            hat = hat .+ MF.estimators(data, floor(Int,sqrt(Tsimu))) ./ Nsimu
            hat_inf = hat_inf .+ MF.mvw_inf(MF.MarkovChainConnectivity(model,θ), excitatory) ./ Nsimu
        end
    
        sum(abs.(hat .- theo)) < 1e-2
        sum(abs.(hat_inf .- theo)) < 1e-4
    end


    Tsimu = Int(1e4)
    N = 1000
    r₊ = rand()
    excitatory = MF.N2excitatory(N, r₊)
    β = rand()
    λ = .3 + .7*rand()
    p = rand()
    model = MarkovChainModel(β*λ,λ,p)

    data = rand(model, excitatory, Tsimu)
    @test begin
        m,v,w = estimators(data)
        m,v,w1 = estimators(data, 1)
        estimators(data, [0,1]) == (m,v,[w,w1])
    end

    @test MF.fit(MarkovChainModel, data, r₊) == MarkovChainModel(MF.Φ(estimators(data)..., r₊)...)

    @testset "inversion of Ψ" begin    
        for (λ,p) in Iterators.product(.1:.1:.9,.1:.1:.9)
            for μ in .1:.1:λ
                for r₊ in [0, .1, .2, .5, .6, .7, .8, .9, 1]
                    @test sum(abs.(MF.Φ(mvw(μ, λ, p, r₊)..., r₊) .- (μ, λ, p))) < 1e-6
                end
            end
        end
    end
end