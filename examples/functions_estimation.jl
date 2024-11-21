using MeanFieldGraph
using ProgressLogging
using DataFrames
using Random
using Statistics
using StatsPlots
using Unicode
using CSV
using TableMetadataTools

## Simulation and computation of the estimators
"""
    estimatorsarray(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int}; id::Int=1, idmax::Int=1)::Tuple{Array{Float64, 4}, Array{Float64, 2}}

Simulate 'Nsimu' times the prescribed model and compute each time the estimators ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` at the observations times ``T`` given by the increasing vector `tvec` (and ``T\\to \\infty``) and the tuning parameters ``\\Delta`` given by `Δvec`.

# Output
Two `Array`s :
1. Its dimension is `6`x`length(Δvec)`x`length(tvec)`x`Nsimu`: the first index correspond (in this order) to the estimators ``\\hat{m}``, ``\\hat{v}``, ``\\hat{w}``, ``\\hat{\\mu}``, ``\\hat{\\lambda}`` and ``\\hat{p}``.
2. Its dimension is `6`x`Nsimu`: as before but it is the limit value when ``T = \\infty``.

# Keywords
They are only used for progress logging.
"""
function estimatorsarray(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int}; id::Int=1, idmax::Int=1)::Tuple{Array{Float64, 4}, Array{Float64, 2}}
    E = Array{Float64}(undef, 6, length(Δvec), length(tvec), Nsimu)
    E_inf = Array{Float64}(undef, 6, Nsimu)
    
    excitatory = MeanFieldGraph.N2excitatory(N, r₊)
    T = tvec[end]

    @progress "estimatorstable, i="*string(id)*" on "*string(idmax) for idsimu in 1:Nsimu
        θ = rand(MeanFieldGraph.ErdosRenyiGraph(N, model.p))
        modelconnec = MarkovChainConnectivity(model,θ)

        m̂, v̂, ŵ = mvw_inf(modelconnec, excitatory)
        μλphat = MeanFieldGraph.Φ(m̂, v̂, ŵ, r₊)
        E_inf[:,idsimu] = [m̂, v̂, ŵ, μλphat...]

        data = rand(modelconnec, excitatory, T)
        for idt in eachindex(tvec)
            tmpdata = data[1:tvec[idt]]
            m̂, v̂, ŵ = estimators(tmpdata, Δvec)
            for idΔ in eachindex(Δvec)
                μλphat = MeanFieldGraph.Φ(m̂, v̂, ŵ[idΔ], r₊)
                E[:,idΔ,idt,idsimu] = [m̂, v̂, ŵ[idΔ], μλphat...]
            end
        end
    end
    return E, E_inf
end

function estimatorstable(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Δ::Int64, Nsimu::Int64}, tvec)::Tuple{DataFrame, DataFrame}
    ## Set default
    N, r₊, β, λ, p, Δ, Nsimu = default_values
    Paramsymbol in (:N, :Δ) ? Typeparam = Int : Typeparam = Float64

    if Paramsymbol == :Δ
        model = MarkovChainModel(β*λ,λ,p)
        E, E_inf = estimatorsarray(model, N, r₊, Nsimu, tvec, Paramvec)
        E = permutedims(E, (1,3,4,2))
        E_inf = reshape(E_inf, (6,Nsimu,1))
    else
        E = Array{Float64}(undef, 6, length(tvec), Nsimu, length(Paramvec))
        E_inf = Array{Float64}(undef, 6, Nsimu, length(Paramvec))
        for idparam in eachindex(Paramvec)
            if Paramsymbol == :N
                N = Paramvec[idparam]
            elseif Paramsymbol == :r₊
                r₊ = Paramvec[idparam]
            elseif Paramsymbol == :β
                β = Paramvec[idparam]
            elseif Paramsymbol == :λ
                λ = Paramvec[idparam]
            elseif Paramsymbol == :p
                p = Paramvec[idparam]
            end
            model = MarkovChainModel(β*λ,λ,p)
            E_new, E_inf_new = estimatorsarray(model, N, r₊, Nsimu, tvec, [Δ]; id=idparam, idmax=length(Paramvec))
            E[:,:,:,idparam] = E_new[:,1,:,:]
            E_inf[:,:,idparam] = E_inf_new
        end
    end

    df = DataFrame(T = Int[], parameter = Typeparam[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[], μ̂ = Float64[], λ̂ = Float64[], p̂ = Float64[])
    complete!(df, E, Paramsymbol, Paramvec, default_values, tvec)
    df_inf = DataFrame(parameter = Typeparam[], m̂ = Float64[], v̂ = Float64[], ŵ = Float64[], μ̂ = Float64[], λ̂ = Float64[], p̂ = Float64[])
    complete!(df_inf, E_inf, Paramsymbol, Paramvec, default_values)
    
    return df, df_inf
end

function metadatacomplete!(df::DataFrame, Paramsymbol::Symbol, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Δ::Int64, Nsimu::Int64})
    N, r₊, β, λ, p, Δ, Nsimu = default_values

    metadata!(df, "Caption", "Values of the estimators")
    Paramsymbol == :N || metadata!(df, "N", N, style=:note)
    Paramsymbol == :r₊ || metadata!(df, "r₊", r₊, style=:note)
    Paramsymbol == :β || metadata!(df, "β", β, style=:note)
    Paramsymbol == :λ || metadata!(df, "λ", λ, style=:note)
    Paramsymbol == :p || metadata!(df, "p", p, style=:note)
    metadata!(df, "Number of simulations", Nsimu, style=:note)
    if "T" in names(df)
        if Paramsymbol != :Δ
            metadata!(df, "Δ", Δ, style=:note)
        end
        colmetadata!(df, :T, "label", "Time horizon used for the estimation", style=:note)
    end
    metadata!(df, "Varying parameter", String(Paramsymbol), style=:note)
    colmetadata!(df, :parameter, "label", "Value of the varying parameter", style=:note)
end

function complete!(df, E, Paramsymbol, Paramvec, default_values, tvec)
    n1, n2, n3, n4 = size(E)
    for idparam in 1:n4
        for idsimu in 1:n3
            for idt in 1:n2
                push!(df, (tvec[idt], Paramvec[idparam], E[:,idt,idsimu,idparam]...))
            end
        end
    end
    metadatacomplete!(df, Paramsymbol, default_values)
end

function complete!(df_inf, E_inf, Paramsymbol, Paramvec, default_values)
    n1, n2, n3 = size(E_inf)
    for idparam in 1:n3
        for idsimu in 1:n2
            push!(df_inf, (Paramvec[idparam], E_inf[:,idsimu,idparam]...))
        end
    end
    metadatacomplete!(df_inf, Paramsymbol, default_values)
end

## Post simulation operations : computation of the errors, means, quantiles, etc.
function estimators2errors(df::DataFrame)::DataFrame
    output = deepcopy(df)
    b = (size(df)[2] == 8)
    paramstring = metadata(df,"Varying parameter")

    paramstring == "β" || (β = metadata(df,"β"))
    paramstring == "λ" || (λ = metadata(df,"λ"))
    paramstring == "p" || (p = metadata(df,"p"))
    paramstring == "r₊" || (r₊ = metadata(df, "r₊"))
    if !(paramstring in ["β", "λ", "p", "r₊"])
        model = MarkovChainModel(β*λ, λ, p)
        m, v, w = mvw(model, r₊)
        output.m̂ = abs.(df.m̂ .- m)
        output.v̂ = abs.(df.v̂ .- v)
        output.ŵ = abs.(df.ŵ .- w)
        output.μ̂ = abs.(df.μ̂ .- model.μ)
        output.λ̂ = abs.(df.λ̂ .- model.λ)
        output.p̂ = abs.(df.p̂ .- model.p)
    else
        for param in unique(df.parameter)
            paramstring == "β" && (β = param)
            paramstring == "λ" && (λ = param)
            paramstring == "p" && (p = param)
            paramstring == "r₊" && (r₊ = param)

            model = MarkovChainModel(β*λ, λ, p)
            m, v, w = mvw(model, r₊)
            selection = (df.parameter .== param)
            output[selection,2+b] = abs.(df[selection,:].m̂ .- m)
            output[selection,3+b] = abs.(df[selection,:].v̂ .- v)
            output[selection,4+b] = abs.(df[selection,:].ŵ .- w)
            output[selection,5+b] = abs.(df[selection,:].μ̂ .- model.μ)
            output[selection,6+b] = abs.(df[selection,:].λ̂ .- model.λ)
            output[selection,7+b] = abs.(df[selection,:].p̂ .- model.p)
        end
    end

    metadata!(output, "Caption", "Values of the absolute difference between the estimator and the theoric value")
    return output
end

function columnquantile(df::DataFrame, q::Real)
    output = empty(df)

    if !("T" in names(df))
        for param in unique(df.parameter)
            selection = (df.parameter .== param)
            quantiles_new = []
            for k in 1:6
                push!(quantiles_new, quantile(df[selection,1+k], q))
            end
            push!(output, (param, quantiles_new...))
        end
    else
        for param in unique(df.parameter)
            for t in unique(df.T) 
                selection = (df.parameter .== param) .* (df.T .== t)
                quantiles_new = []
                for k in 1:6
                    push!(quantiles_new, quantile(df[selection,2+k], q))
                end
                push!(output, (t, param, quantiles_new...))
            end
        end
    end    

    metadata!(output, "Caption", string(100*q)*"% quantile of the absolute difference between the estimator and the theoric value")
    return output
end

## Plot functions
function plotestimators(df::DataFrame, df_inf::DataFrame; targetstring::String="all")
    if targetstring == "all"
        for ts in ["m̂","v̂","ŵ","μ̂","λ̂","p̂"]
            plotestimators(df, df_inf; targetstring = ts)
        end
    else
        targetstring = Unicode.normalize(targetstring)
        paramvec = unique(df.parameter)
        paramstring = metadata(df, "Varying parameter")
    
        plot(df.T, df[:,targetstring], group=df.parameter, color=indexin(df.parameter,paramvec))
        scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring], group=df_inf.parameter, color=indexin(df_inf.parameter,paramvec), marker = :hline, label=false, markerstrokewidth = 5)
        xlabel!("T")
        ylims!(0, min(maximum(df[:,targetstring]),.3))
        ylabel!("Absolute error")
        display(title!("Estimation error for "*targetstring*" as "*paramstring*" varies, Δ="*string(metadata(df,"Δ"))))
    end
end

function plotestimators(dfs::Tuple{DataFrame,DataFrame,DataFrame}, dfs_inf::Tuple{DataFrame,DataFrame,DataFrame}; targetstring::String="all")
    if targetstring == "all"
        for ts in ["m̂","v̂","ŵ","μ̂","λ̂","p̂"]
            plotestimators(dfs, dfs_inf; targetstring = ts)
        end
    else
        dflow, df, dfup = dfs
        dflow_inf, df_inf, dfup_inf = dfs_inf
        targetstring = Unicode.normalize(targetstring)
        paramvec = unique(df.parameter)
        paramstring = metadata(df, "Varying parameter")
    
        Tmax = maximum(df.T)
        dflowup_inf = empty(df)
        for param in paramvec
            push!(dflowup_inf, (Tmax, param, Array(dflow_inf[dflow_inf.parameter .== param, Not(1)])...))
            push!(dflowup_inf, (Tmax, param, Array(dfup_inf[dfup_inf.parameter .== param, Not(1)])...))
        end
        plot(df.T, df[:,targetstring], group=df.parameter, color=indexin(df.parameter,paramvec); ribbon = (dflow[:,targetstring], dfup[:,targetstring]))
        scatter!(maximum(df.T)*ones(size(df_inf)[1]), df_inf[:,targetstring], group=df_inf.parameter, color=indexin(df_inf.parameter,paramvec), marker = :hline, label=false, markerstrokewidth = 2)
        plot!(dflowup_inf.T, dflowup_inf[:,targetstring], group=dflowup_inf.parameter, color=indexin(dflowup_inf.parameter,paramvec), label=false)
        xlabel!("T")
        ylims!(0, min(maximum(df[:,targetstring]),.3))
        ylabel!("Absolute error")
        display(title!("Estimation error for "*targetstring*" as "*paramstring*" varies, Δ="*string(metadata(df,"Δ"))))
    end
end

## Save and Load functions
function simulationandsave(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Δ::Int64, Nsimu::Int64}, tvec)
    df, df_inf = estimatorstable(Paramsymbol, Paramvec, default_values, tvec)
    paramstring = string(Paramsymbol)

    CSV.write("data/estimators_vary_"*paramstring*"_delta_"*string(default_values.Δ)*".csv", df)
    open("data/estimators_vary_"*paramstring*"_delta_"*string(default_values.Δ)*".toml", "w") do io
        print(io, meta2toml(df))
    end

    CSV.write("data/estimators_vary_"*paramstring*"_delta_"*string(default_values.Δ)*"_inf.csv", df_inf)
    open("data/estimators_vary_"*paramstring*"_delta_"*string(default_values.Δ)*"_inf.toml", "w") do io
        print(io, meta2toml(df_inf))
    end
end

function estimatorsload(filename::String)
    df = CSV.read(filename*".csv", DataFrame)
    open(filename*".toml") do io
        toml2meta!(df, io)
    end
    df_inf = CSV.read(filename*"_inf.csv", DataFrame)
    open(filename*"_inf.toml") do io
        toml2meta!(df_inf, io)
    end

    return df, df_inf
end