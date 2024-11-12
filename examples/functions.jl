using MeanFieldGraph
using ProgressLogging
using DataFrames
using Random
using Statistics
using StatsPlots

"""
    estimatorstable(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int})::Tuple{DataFrame, DataFrame}

Simulate 'Nsimu' times the prescribed model and compute each time the estimators ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` at the observations times ``T`` given by `tvec` (and ``T\\to \\infty``) and the tuning parameters ``\\Delta`` given by `Δvec`.

# Output
Two `DataFrame` :
1. For each simulation correspond `length(tvec)` rows with the computed values of ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` at the prescribed times ``T``.
2. For each simulation corresponds one row with the computed values of ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` with ``T = \\infty``.
"""
function estimatorstable(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int})::Tuple{DataFrame, DataFrame}
    T = maximum(tvec)
    df = DataFrame(idsimu = Int[], T = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Vector{Float64}[])
    df_inf = DataFrame(idsimu = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[])
    
    excitatory = MeanFieldGraph.N2excitatory(N, r₊)
    
    @progress "simulation" for idsimu in 1:Nsimu
        θ = rand(MeanFieldGraph.ErdosRenyiGraph(N, model.p))
        modelconnec = MarkovChainConnectivity(model,θ)
        push!(df_inf, (idsimu, mvw_inf(modelconnec, excitatory)...))

        data = rand(modelconnec, excitatory, T)
        for t in tvec
            tmpdata = data[1:t]
            push!(df, (idsimu, t, estimators(tmpdata, Δvec)...))
        end
    end
    return df, df_inf
end

function extractw_computeμλp(df::DataFrame, column::Int)::DataFrame
    if "T" in names(df)
        output = DataFrame(idsimu = Int[], T = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[], μ̂ = Float64[], λ̂ = Float64[], p̂ = Float64[])
    else
        output = DataFrame(idsimu = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[], μ̂ = Float64[], λ̂ = Float64[], p̂ = Float64[])
    end

    for id in 1:nrow(df)
        m̂ = df.m̂[id]
        v̂ = df.v̂[id]
        ŵ = df.ŵ[id][column]
        μλphat = MeanFieldGraph.Φ(m̂, v̂, ŵ, r₊)
        if "T" in names(df)
            push!(output, (df.idsimu[id], df.T[id], m̂, v̂, ŵ, μλphat...))
        else
            push!(output, (df.idsimu[id], m̂, v̂, ŵ, μλphat...))
        end
    end

    return output
end

function estimators2errors!(df::DataFrame, model::MarkovChainModel, r₊::Real)
    m, v, w = mvw(model, r₊)
    df.m̂ = abs.(df.m̂ .- m)
    df.v̂ = abs.(df.v̂ .- v)
    df.ŵ = abs.(df.ŵ .- w)
    df.μ̂ = abs.(df.μ̂ .- model.μ)
    df.λ̂ = abs.(df.λ̂ .- model.λ)
    df.p̂ = abs.(df.p̂ .- model.p)
end

function columnquantile(df::DataFrame, q::Real)
    if !("T" in names(df))
        return DataFrame(m̂ = quantile(df.m̂, q),
                            v̂ = quantile(df.v̂, q),
                            ŵ = quantile(df.ŵ, q),
                            μ̂ = quantile(df.μ̂, q),
                            λ̂ = quantile(df.λ̂, q),
                            p̂ = quantile(df.p̂, q))
    end
    
    output = DataFrame(T = Int[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[], μ̂ = Float64[], λ̂ = Float64[], p̂ = Float64[])

    for T in unique(df.T)
        push!(output, (T,    quantile(df.m̂[df.T .== T], q),
                            quantile(df.v̂[df.T .== T], q),
                            quantile(df.ŵ[df.T .== T], q),
                            quantile(df.μ̂[df.T .== T], q),
                            quantile(df.λ̂[df.T .== T], q),
                            quantile(df.p̂[df.T .== T], q)))
    end

    return output
end