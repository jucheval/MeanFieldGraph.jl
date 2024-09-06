using DiscreteHawkes
using DataFrames
using Distributions
using ProgressLogging
using StatsPlots

# model
N = 1000
r₊ = .7
β = .5
λ = .05
p = .95
μ = β*λ
model = Model(N,r₊,β,λ,p)

# parameters of the method
T = 5000
Δ = 70 # ≈ sqrt(T)
Nsimu = 100

# estimation
begin 
    targets = [μ,λ,p]
    df_estim = DataFrame(simulationid = Int[], estimatorname = String[], estimatorid = Int[], estimatorvalue = Float64[])
    parameternames = ["μ̂", "λ̂", "p̂"]
    @progress "simulation and estimation" for simulationid in 1:Nsimu
        data = rand(model, T)
        estimators = DiscreteHawkes.fit(data)
        if typeof(estimators) != Tuple{Float64, Float64, Float64}
            estimators = estimators[1][2:4]
        end
        for index in 1:3
            push!(df_estim, (simulationid, parameternames[index], index, estimators[index]))
        end
    end
end

# plot
begin 
    p1 = violin(df_estim.estimatorvalue[df_estim.estimatorid .== 1], label="")
    scatter!(p1, [(1, targets[1])], color=:tomato, label="μ")
    p2 = violin(df_estim.estimatorvalue[df_estim.estimatorid .== 2], label="")
    scatter!(p2, [(1, targets[2])], color=:tomato, label="λ")
    p3 = violin(df_estim.estimatorvalue[df_estim.estimatorid .== 3], label="")
    scatter!(p3, [(1, targets[3])], color=:tomato, label="p")
    plot(p1, p2, p3, layout = (1,3))
end


