include("functions_classification.jl")
Random.seed!(1)

## Default values
r₊ = .5
β = .5
λ = .5
p = .5

## Specific values
Nvec = 50:50:600
T = Int(1e5)
length_tvec = 20
tmin = 1000
tvec = floor.(Int,collect(range(tmin, T, length_tvec)))
Nsimu = 100


model = MarkovChainModel(β*λ,λ,p)
begin
    X = zeros(length(tvec), length(Nvec))
    for idN in eachindex(Nvec)
        excitatory = MeanFieldGraph.N2excitatory(Nvec[idN], r₊)
        @progress "table for color plot, i="*string(idN)*" on "*string(length(Nvec)) for idsimu in 1:Nsimu
            data = rand(model, excitatory, T)
            for idt in eachindex(tvec)
                tmpdata = data[1:tvec[idt]]
                naive, kmeans = classification(tmpdata)
                X[idt, idN] += (kmeans == excitatory)
            end
        end
    end
    X = X/Nsimu
end

using GLMakie
fig = Figure()
ax = Axis(fig[1, 1], xlabel = "T", ylabel = "N")
hm = Makie.heatmap!(tvec, Nvec, X, interpolate = true)
