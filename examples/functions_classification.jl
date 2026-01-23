using MeanFieldGraph
using Clustering
using KrylovKit
using ProgressLogging
using DataFrames
using Random
using Statistics
using StatsPlots
using Distributions
using CSV
using TableMetadataTools
using LaTeXStrings

# FIXME : since σ̂sp is estimated up to +/- 1 factor, we use σ̂ag to decide which cluster is excitatory and which is inhibitory

#region Simulation and creation of the dataset
"""
    classiferrortable(Paramsymbol::Symbol, Paramvec, default_values::@NamedTuple{N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64}, tvec)::DataFrame

For each value of the parameter `Paramsymbol` prescribed via `Paramvec`, make `default_values.Nsimu` simulations of the model, compute the misclassification proportion at the observations times ``T`` given by `tvec`. The fixed parameters are specified by the named tuple `default_values`.

# Output
One `DataFrame`: each parameter value correspond to `default_values.Nsimu * length(tvec)` rows. The last 10 columns give the proportion of misclassification for each coupling of method (aggregated or spectral) and clustering (kmeans or hierarchical with complete, single, average, ward linkage).
"""
function classiferrortable(
    Paramsymbol::Symbol,
    Paramvec,
    default_values::@NamedTuple{
        N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64
    },
    tvec;
    multi_thread::Bool=false,
)::DataFrame
    ## Set default
    N, r₊, β, λ, p, Nsimu = default_values

    # FIXME A = Array{Float64}(undef, 10, length(tvec), Nsimu, length(Paramvec))
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
        model = MarkovChainModel(β * λ, λ, p)
        excitatory = MeanFieldGraph.N2excitatory(N, r₊)
        if multi_thread
            println(
                "classiferrortables, i=" *
                string(idparam) *
                " on " *
                string(length(Paramvec)),
            )
            Threads.@threads for idsimu in 1:Nsimu
                data = rand(model, excitatory, T)
                for idt in eachindex(tvec)
                    tmpdata = data[1:tvec[idt]]
                    # FIXME A[:, idt, idsimu, idparam] = misclassificationrates(tmpdata, excitatory)
                    σ̂ag = MeanFieldGraph.covariance_vector(tmpdata)
                    A[idt, idsimu, idparam] = misclassificationrate(
                        σ̂ag, excitatory, :kmeans
                    )
                end
            end
        else
            @progress "classiferrortables, i=" *
                string(idparam) *
                " on " *
                string(length(Paramvec)) for idsimu in 1:Nsimu
                data = rand(model, excitatory, T)
                for idt in eachindex(tvec)
                    tmpdata = data[1:tvec[idt]]
                    A[:, idt, idsimu, idparam] = misclassificationrates(tmpdata, excitatory)
                end
            end
        end
    end

    Paramsymbol == :N ? Typeparam = Int : Typeparam = Float64
    df = DataFrame(;
        parameter=Typeparam[],
        T=Int[],
        miscl_rate_ag_kmeans=Float64[],
        # miscl_rate_ag_hccomp=Float64[], FIXME
        # miscl_rate_ag_hcsing=Float64[],
        # miscl_rate_ag_hcavg=Float64[],
        # miscl_rate_ag_hcward=Float64[],
        # miscl_rate_sp_kmeans=Float64[],
        # miscl_rate_sp_hccomp=Float64[],
        # miscl_rate_sp_hcsing=Float64[],
        # miscl_rate_sp_hcavg=Float64[],
        # miscl_rate_sp_hcward=Float64[],
    )
    # _, n1, n2, n3 = size(A) FIXME
    n1, n2, n3 = size(A)
    for idparam in 1:n3
        for idsimu in 1:n2
            for idt in 1:n1
                push!(df, (
                    Paramvec[idparam],
                    tvec[idt],
                    # A[:, idt, idsimu, idparam]...FIXME
                    A[idt, idsimu, idparam],
                ))
            end
        end
    end
    metadatacomplete!(df, Paramsymbol, default_values)
    return df
end

#region Helper functions used in classiferrortable
## fill the metadata of the DataFrame
function metadatacomplete!(
    df::DataFrame,
    Paramsymbol::Symbol,
    default_values::@NamedTuple{
        N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64
    }
)
    N, r₊, β, λ, p, Nsimu = default_values

    metadata!(df, "Caption", "Misclassification rates")
    Paramsymbol == :N || metadata!(df, "N", N; style=:note)
    Paramsymbol == :r₊ || metadata!(df, "r₊", r₊; style=:note)
    Paramsymbol == :β || metadata!(df, "β", β; style=:note)
    Paramsymbol == :λ || metadata!(df, "λ", λ; style=:note)
    Paramsymbol == :p || metadata!(df, "p", p; style=:note)
    metadata!(df, "Number of simulations", Nsimu; style=:note)
    metadata!(df, "Varying parameter", String(Paramsymbol); style=:note)

    colmetadata!(df, :T, "label", "Time horizon used for the estimation"; style=:note)
    methods_str = ["aggregated", "spectral"]
    methods_symb = [:ag, :sp]
    clusterings_str = [
        "kmeans",
        "hierarchical (complete)",
        "hierarchical (single)",
        "hierarchical (average)",
        "hierarchical (ward)",
    ]
    clusterings_symb = [:kmeans, :hccomp, :hcsing, :hcavg, :hcward]
    # for i in 1:2 FIXME
    #     for j in 1:5
    #         colmetadata!(
    #             df,
    #             Symbol(:miscl_rate_, methods_symb[i], :_, clusterings_symb[j]),
    #             "label",
    #             "Misclassification rate for the " *
    #             methods_str[i] *
    #             " method with " *
    #             clusterings_str[j] *
    #             " clustering";
    #             style=:note,
    #         )
    #     end
    # end
    colmetadata!(
        df,
        :miscl_rate_ag_kmeans,
        "label",
        "Misclassification rate for the aggregated method with kmeans clustering";
        style=:note,
    )
    return colmetadata!(
        df, :parameter, "label", "Value of the varying parameter"; style=:note
    )
end

## Compute misclassification rate for all methods
function misclassificationrates(
    data::DiscreteTimeData, excitatory::Vector{Bool}
)::Vector{Float64}
    output = Float64[]

    # aggregated classifications
    σ̂ag = MeanFieldGraph.covariance_vector(data)
    for c in [:kmeans, :complete, :single, :average, :ward]
        push!(output, misclassificationrate(σ̂ag, excitatory, c))
    end

    # spectral classifications
    Σ̂ = MeanFieldGraph.covariance_matrix(data)
    _, vecs = eigsolve(transpose(Σ̂) * Σ̂)
    σ̂sp = vecs[1]
    # FIXME : how to chose the sign of σ̂sp ?
    if mapreduce(abs, +, σ̂ag - σ̂sp) > mapreduce(abs, +, σ̂ag + σ̂sp)
        σ̂sp *= -1
    end
    for c in [:kmeans, :complete, :single, :average, :ward]
        push!(output, misclassificationrate(σ̂sp, excitatory, c))
    end

    return output
end

## Compute misclassification rate for each individual methods
function misclassificationrate(
    σ::Vector{Float64}, excitatory::Vector{Bool}, clustering::Symbol
)::Float64
    if clustering == :kmeans
        classif = MeanFieldGraph.cluster2bool(
            kmeans(transpose(σ), 2; init=[argmin(σ), argmax(σ)])
        )
    else
        distances = [abs(σ[i] - σ[j]) for i in eachindex(σ), j in eachindex(σ)]
        ct = cutree(hclust(distances; linkage=clustering); k=2)
        id_excitatory = ct[argmax(σ)]
        classif = ct .== id_excitatory
    end

    return mean(abs.(classif - excitatory))
end
#endregion
#endregion

# TODO: restart here. Leave the dataframe as it is. Make a function that makes the dataframe tidy afterwards

# TODO : convert into Tidier syntax (maybe simpler to move it to the "call" files)
#region Post simulation operations
## Compute the mean misclassification rate and probability of exact recovery
function proportion2errors(df::DataFrame)
    output = empty(df)
    select!(output, Not(3))
    output[!, :mean_prop] .= 0.0
    output[!, :mean_prop_std] .= 0.0
    output[!, :proba_exact_recovery] .= 0.0
    output[!, :proba_exact_recovery_std] .= 0.0

    for param in unique(df.parameter)
        for t in unique(df.T)
            selection = (df.parameter .== param) .* (df.T .== t)
            push!(
                output,
                (
                    param,
                    t,
                    mean(df.prop_errors[selection]),
                    std(df.prop_errors[selection]),
                    mean(df.prop_errors[selection] .== 0),
                    std(df.prop_errors[selection] .== 0),
                ),
            )
        end
    end
    metadata!(
        output,
        "Caption",
        "Mean proportion of misclassification and probability of exact recovery",
    )
    return output
end

### Compute the mean misclassification rate and probability of exact recovery + bands for a wide data set
function mmr_per(df_wide)
    # lengthen data set and select columns
    df = @chain df_wide begin
        @rename(n = parameter, t = T)
        @pivot_longer(-[n, t], values_to = "misclassification_rate")
        @separate(variable, (a, b, method, clustering), "_")
        @select(-[a, b])
    end

    # compute mean and standard error for misclassification rate and exact recovery
    df_mean = @chain df begin
        @group_by([n, t, method, clustering])
        @summarize(
            mr_mean = mean(misclassification_rate),
            mr_std = std(misclassification_rate),
            er_mean = mean(misclassification_rate == 0),
            er_std = std(misclassification_rate == 0)
        )
        @ungroup
    end

    # compute lower and upper bounds for bands
    factor = quantile(Normal(), 0.975) / sqrt(metadata(df_mean, "Number of simulations"))
    df_mean_bands = transform(
        df_mean,
        [:mr_mean, :mr_std] => ((x, y) -> x - y) => :mr_lower,
        [:mr_mean, :mr_std] => ((x, y) -> x + y) => :mr_upper,
        [:er_mean, :er_std] => ((x, y) -> x - factor * y) => :er_lower,
        [:er_mean, :er_std] => ((x, y) -> x + factor * y) => :er_upper,
    )

    return df_mean_bands
end
#endregion

# TODO: Convert to AlgebraofGraphics plots ?? 
#region Plot functions
function plotclassification(df::DataFrame)
    paramvec = unique(df.parameter)
    paramstring = metadata(df, "Varying parameter")
    Nsimu = metadata(df, "Number of simulations")
    col = indexin(df.parameter, paramvec)
    q = quantile(Normal(), 0.975)

    paramstring == "r₊" && (paramstring = "r_+") # modified before passing to latexstring

    plot(df.T, df.mean_prop; group=df.parameter, color=col, ribbon=df.mean_prop_std)
    xlabel!(L"T")
    ylims!(0, 1)
    title!(
        "Misclassification and exact recovery as " * latexstring(paramstring) * " varies"
    )
    return display(
        plot!(
            df.T,
            df.proba_exact_recovery;
            group=df.parameter,
            color=col,
            label=false,
            linestyle=:dash,
            ribbon=q * df.proba_exact_recovery_std / sqrt(Nsimu),
        ),
    )
end
#endregion

#region Save and Load functions
function simulationandsave(
    Paramsymbol::Symbol,
    Paramvec,
    default_values::@NamedTuple{
        N::Int64, r₊::Float64, β::Float64, λ::Float64, p::Float64, Nsimu::Int64
    },
    tvec,
)
    df = classiferrortable(Paramsymbol, Paramvec, default_values, tvec)
    paramstring = string(Paramsymbol)

    CSV.write("data/CO24/classification_vary_" * paramstring * ".csv", df)
    open("data/CO24/classification_vary_" * paramstring * ".toml", "w") do io
        print(io, meta2toml(df))
    end
end

function estimatorsload(filename::String)
    df = CSV.read(filename * ".csv", DataFrame)
    open(filename * ".toml") do io
        toml2meta!(df, io)
    end

    if !(isfile(filename * "_inf.csv"))
        return df
    end
    df_inf = CSV.read(filename * "_inf.csv", DataFrame)
    open(filename * "_inf.toml") do io
        toml2meta!(df_inf, io)
    end

    return df, df_inf
end
#endregion