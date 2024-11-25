using MeanFieldGraph
using ProgressLogging
using DataFrames
using Random
using Statistics
using StatsPlots
using Distributions
using CSV
using TableMetadataTools
using LaTeXStrings

## Simulation and computation of the estimators
"""
    classification_errors(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int}, Δvec::Vector{Int})::Tuple{DataFrame, DataFrame}

Simulate 'Nsimu' times the prescribed model and compute each time the proportion of classification errors at the observations times ``T`` given by `tvec`.

# Output
One `DataFrame`: each simulation correspond to `length(tvec)` rows. The three columns give the index of the simulation, the observation time and the proportion of classification errors.
"""
function classification_errors(model::MarkovChainModel, N::Int, r₊::Float64, Nsimu::Int, tvec::Vector{Int})::DataFrame
    T = maximum(tvec)
    df = DataFrame(idsimu = Int[], T = Int[], prop_errors = Float64[])
    
    excitatory = MeanFieldGraph.N2excitatory(N, r₊)
    
    @progress "classification_errors" for idsimu in 1:Nsimu
        data = rand(model, excitatory, T)
        for t in tvec
            tmpdata = data[1:t]
            partition = classification(tmpdata)
            push!(df, (idsimu, t, mean(abs.(partition - excitatory))))
        end
    end
    return df
end

"""
    classiferrortable(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64}, tvec)::DataFrame

For each value of the parameter `Paramsymbol` prescribed via `Paramvec`, make `default_values.Nsimu` simulations of the mmodel, compute the classification at the observations times ``T`` given by `tvec` and the misclassification proportion. The fixed parameters are specified by the named tuple `default_values`.

# Output
One `DataFrame`: each parameter value correspond to `default_values.Nsimu * length(tvec)` rows. The last column gives the proportion of misclassification.
"""
function classiferrortable(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64}, tvec)::DataFrame
    ## Set default
    N, r₊, β, λ, p, Nsimu = default_values

    A = Array{Float64}(undef, length(tvec), Nsimu, length(Paramvec))

    for idparam in eachindex(Paramvec)
        ## Set varying parameter
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
        excitatory = MeanFieldGraph.N2excitatory(N, r₊)
        @progress "classiferrortables, i="*string(idparam)*" on "*string(length(Paramvec)) for idsimu in 1:Nsimu
        data = rand(model, excitatory, T)
            for idt in eachindex(tvec)
                tmpdata = data[1:tvec[idt]]
                A[idt, idsimu, idparam] = mean(abs.(classification(tmpdata) - excitatory))
            end
        end
    end

    Paramsymbol == :N ? Typeparam = Int : Typeparam = Float64
    df = DataFrame(parameter = Typeparam[], T = Int[], prop_errors = Float64[])
    n1, n2, n3 = size(A)
    for idparam in 1:n3
        for idsimu in 1:n2 
            for idt in 1:n1
                push!(df, (Paramvec[idparam], tvec[idt], A[idt,idsimu,idparam]))
            end
        end
    end
    metadatacomplete!(df, Paramsymbol, default_values)
    return df
end

function metadatacomplete!(df::DataFrame, Paramsymbol::Symbol, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64})
    N, r₊, β, λ, p, Nsimu = default_values

    metadata!(df, "Caption", "Proportion of misclassification")
    Paramsymbol == :N || metadata!(df, "N", N, style=:note)
    Paramsymbol == :r₊ || metadata!(df, "r₊", r₊, style=:note)
    Paramsymbol == :β || metadata!(df, "β", β, style=:note)
    Paramsymbol == :λ || metadata!(df, "λ", λ, style=:note)
    Paramsymbol == :p || metadata!(df, "p", p, style=:note)
    metadata!(df, "Number of simulations", Nsimu, style=:note)
    colmetadata!(df, :T, "label", "Time horizon used for the estimation", style=:note)
    metadata!(df, "Varying parameter", String(Paramsymbol), style=:note)
    colmetadata!(df, :parameter, "label", "Value of the varying parameter", style=:note)
end

## Post simulation operations : computation of the errors
function proportion2errors(df::DataFrame)
    output = empty(df)
    select!(output, Not(3))
    output[!,:mean_prop] .= 0.
    output[!,:mean_prop_std] .= 0.
    output[!,:proba_exact_recovery] .= 0.
    output[!,:proba_exact_recovery_std] .= 0.

    for param in unique(df.parameter)
        for t in unique(df.T) 
            selection = (df.parameter .== param) .* (df.T .== t)
            push!(output, (param, t, mean(df.prop_errors[selection]), std(df.prop_errors[selection]), mean(df.prop_errors[selection] .== 0), std(df.prop_errors[selection] .== 0)))
        end
    end
    metadata!(output, "Caption", "Mean proportion of misclassification and probability of exact recovery")
    return output
end

## Plot functions
function plotclassification(df::DataFrame)
    paramvec = unique(df.parameter)
    paramstring = metadata(df, "Varying parameter")
    Nsimu = metadata(df, "Number of simulations")
    col = indexin(df.parameter,paramvec)
    q = quantile(Normal(), 0.975)
    
    paramstring == "r₊" && (paramstring = "r_+") # modified before passing to latexstring

    plot(df.T, df.mean_prop, group=df.parameter, color=col; ribbon = df.mean_prop_std)
    xlabel!(L"T")
    ylims!(0,1)
    title!("Misclassification and exact recovery as "*latexstring(paramstring)*" varies")
    display(plot!(df.T, df.proba_exact_recovery, group=df.parameter, color=col, label=false, linestyle = :dash; ribbon = q*df.proba_exact_recovery_std/sqrt(Nsimu)))
end

## Save and Load functions
function simulationandsave(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64}, tvec)
    df = classiferrortable(Paramsymbol, Paramvec, default_values, tvec)
    paramstring = string(Paramsymbol)

    CSV.write("data/CO24/classification_vary_"*paramstring*".csv", df)
    open("data/CO24/classification_vary_"*paramstring*".toml", "w") do io
        print(io, meta2toml(df))
    end
end

function estimatorsload(filename::String)
    df = CSV.read(filename*".csv", DataFrame)
    open(filename*".toml") do io
        toml2meta!(df, io)
    end

    if !(isfile(filename*"_inf.csv"))
        return df
    end
    df_inf = CSV.read(filename*"_inf.csv", DataFrame)
    open(filename*"_inf.toml") do io
        toml2meta!(df_inf, io)
    end

    return df, df_inf
end