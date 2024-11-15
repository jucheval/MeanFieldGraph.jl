using MeanFieldGraph
using ProgressLogging
using DataFrames
using Random
using Statistics
using StatsPlots

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

function classiferrorarray(Paramsymbol::Symbol, Paramvec, default_values, Nsimu, tvec)::Array{Float64}
    output = Array{Float64}(undef, length(tvec), Nsimu, length(Paramvec), 2)

    ## Set default
    N, r₊, β, λ, p = default_values
    
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
                naive, kmeans = classification(tmpdata)
                output[idt, idsimu, idparam, 1] = mean(abs.(naive - excitatory))
                output[idt, idsimu, idparam, 2] = mean(abs.(kmeans - excitatory))
            end
        end
    end
    return output
end

"""
    classiferrortables(Paramsymbol::Symbol, Paramvec, default_values, Nsimu, tvec)::DataFrame

For each value of the parameter `Paramsymbol` prescribed via `Paramvec`, estimate the mean proportion of classification errors and the probability of exact recovery over 'Nsimu' simulations at the observations times ``T`` given by `tvec`.

# Output
One `DataFrame`: each parameter value correspond to `length(tvec)` rows. The last two columns give the mean proportion of classification errors and the probability of exact recovery.
"""
function classiferrortables(Paramsymbol::Symbol, Paramvec, default_values, Nsimu, tvec)::DataFrame
    A = classiferrorarray(Paramsymbol::Symbol, Paramvec, default_values, Nsimu, tvec)

    Paramsymbol == :N ? Typeparam = Int : Typeparam = Float64
    df = DataFrame(color = Int[], parameter = Typeparam[], T = Int[], 
                    prop_errors_naive = Float64[], exact_recovery_naive = Float64[],
                    prop_errors_kmeans = Float64[], exact_recovery_kmeans = Float64[])
    n1, n2, n3, n4 = size(A)
    for idparam in 1:n3
        for idt in 1:n1
            push!(df, (idparam, Paramvec[idparam], tvec[idt], 
                        mean(A[idt,:,idparam,1]), mean(A[idt,:,idparam,1] .== 0),
                        mean(A[idt,:,idparam,2]), mean(A[idt,:,idparam,2] .== 0)))
        end
    end
    metadata!(df, "Varying parameter", String(Paramsymbol))
    return df
end

function plotclassification(df::DataFrame; type::String="none")
    if type == "none"
        for type in ["naive", "kmeans", "kmedoids"]
            plotclassification(df, type)
        end
    else
        plotclassification(df, type)
    end
end

function plotclassification(df::DataFrame, type::String)
    paramstring = metadata(df, "Varying parameter")
    pestring = "prop_errors_"*type
    erstring = "exact_recovery_"*type
    plot(df.T, df[!,pestring], group=df.parameter, color=df.color)
    xlabel!("T")
    title!("Classification error as "*paramstring*" varies - "*type)
    display(plot!(df.T, df[!,erstring], group=df.parameter, color=df.color, label=false, linestyle = :dash))
end