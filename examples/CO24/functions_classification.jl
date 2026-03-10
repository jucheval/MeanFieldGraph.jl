using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using MeanFieldGraph
using Clustering
using KrylovKit

using ProgressLogging

using CSV

using Distributions
using Random
using Statistics

using DataFrames
using DataFramesMeta
using TableMetadataTools

using CairoMakie
using AlgebraOfGraphics
using LaTeXStrings

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
    methods::Vector{Symbol}=[:ag, :sp],
    clusterings::Vector{Symbol}=[:kmeans, :hccomp, :hcsing, :hcavg, :hcward],
)::DataFrame
    ## Set default
    N, r₊, β, λ, p, Nsimu = default_values

    A = Array{Float64}(
        undef, length(methods) * length(clusterings), length(tvec), Nsimu, length(Paramvec)
    )

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
            println("Multi-thread activated on $(Threads.nthreads()) threads")
            println(
                "classiferrortables, i=" *
                string(idparam) *
                " on " *
                string(length(Paramvec)) *
                ". " *
                string(Paramsymbol) *
                " = " *
                string(Paramvec[idparam]),
            )
            Threads.@threads for idsimu in 1:Nsimu
                data = rand(model, excitatory, T)
                for idt in eachindex(tvec)
                    tmpdata = data[1:tvec[idt]]
                    A[:, idt, idsimu, idparam] = misclassificationrates(
                        tmpdata, excitatory; methods=methods, clusterings=clusterings
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
                    A[:, idt, idsimu, idparam] = misclassificationrates(
                        tmpdata, excitatory; methods=methods, clusterings=clusterings
                    )
                end
            end
        end
    end

    Paramsymbol == :N ? Typeparam = Int : Typeparam = Float64
    df = DataFrame(; parameter=Typeparam[], T=Int[])
    column_names = [
        Symbol(:miscl_rate_, m, :_, c) => Float64[] for c in clusterings, m in methods
    ]
    insertcols!(df, column_names...)
    _, n1, n2, n3 = size(A)
    for idparam in 1:n3
        for idsimu in 1:n2
            for idt in 1:n1
                push!(df, (Paramvec[idparam], tvec[idt], A[:, idt, idsimu, idparam]...))
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
    methods_fullstr = Dict("ag" => "aggregated", "sp" => "spectral")
    clusterings_fullstr = Dict(
        "kmeans" => "kmeans",
        "hccomp" => "hierarchical (complete)",
        "hcsing" => "hierarchical (single)",
        "hcavg" => "hierarchical (average)",
        "hcward" => "hierarchical (ward)",
    )
    for (i, col_str) in enumerate(names(df))
        col_str in ("parameter", "T") && continue
        _, _, method, clustering = split(col_str, "_")
        colmetadata!(
            df,
            i,
            "label",
            "Misclassification rate for the " *
            methods_fullstr[method] *
            " method with " *
            clusterings_fullstr[clustering] *
            " clustering";
            style=:note,
        )
    end
    return colmetadata!(
        df, :parameter, "label", "Value of the varying parameter"; style=:note
    )
end

## Compute misclassification rate for all methods
function misclassificationrates(
    data::DiscreteTimeData,
    excitatory::Vector{Bool};
    methods::Vector{Symbol}=[:ag, :sp],
    clusterings::Vector{Symbol}=[:kmeans, :hccomp, :hcsing, :hcavg, :hcward],
)::Vector{Float64}
    output = Float64[]

    σ̂ag = MeanFieldGraph.covariance_vector(data)

    for method in methods
        if method == :ag
            # aggregated classifications
            σ̂ = σ̂ag
        else
            # spectral classifications
            Σ̂ = MeanFieldGraph.covariance_matrix(data)
            _, vecs = eigsolve(transpose(Σ̂) * Σ̂)
            σ̂ = vecs[1]
            # FIXME : how to chose the sign of σ̂sp ?
            if mapreduce(abs, +, σ̂ag - σ̂) > mapreduce(abs, +, σ̂ag + σ̂)
                σ̂ *= -1
            end
        end
        for c in clusterings
            push!(output, misclassificationrate(σ̂, excitatory, c))
        end
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
        linkages = Dict(
            :hccomp => :complete, :hcsing => :single, :hcavg => :average, :hcward => :ward
        )
        distances = [abs(σ[i] - σ[j]) for i in eachindex(σ), j in eachindex(σ)]
        ct = cutree(hclust(distances; linkage=linkages[clustering]); k=2)
        id_excitatory = ct[argmax(σ)]
        classif = ct .== id_excitatory
    end

    return mean(abs.(classif - excitatory))
end
#endregion
#endregion

#region Post simulation operations
### Compute the mean misclassification rate and probability of exact recovery + bands for a wide data set
function mmr_per(df_wide)
    output = @chain df_wide begin
        ### reshape from wide to long
        stack(Not([:parameter, :T]); value_name=:misclassification_rate)
        ### extract method and clustering from the column names
        @rtransform! _ @astable begin
            (a, b, c, d) = split(:variable, "_")
            :method = c
            :clustering = d
        end
        ### remove the column with the original column names
        @select!(Not(:variable))
        ### compute mean and standard error for misclassification rate and exact recovery
        @groupby([:parameter, :T, :method, :clustering])
        @combine begin
            :mr_mean = mean(:misclassification_rate)
            :mr_std = std(:misclassification_rate)
            :er_mean = mean(:misclassification_rate .== 0)
            :er_std = std(:misclassification_rate .== 0)
        end
        ### compute lower and upper bounds for bands
        @transform! _ @astable begin
            # multiplicative factor for the confidence interval of the exact recovery rate
            a = quantile(Normal(), 0.975) / sqrt(metadata(_, "Number of simulations"))
            :mr_lower = :mr_mean - :mr_std
            :mr_upper = :mr_mean + :mr_std
            :er_lower = :er_mean - a * :er_std
            :er_upper = :er_mean + a * :er_std
        end
    end

    return output
end
#endregion

#region Plot functions
function plotclassification(df::DataFrame)
    paramstring = metadata(df, "Varying parameter")
    paramstring == "r₊" && (paramstring = "r_+") # modified before passing to latexstring

    @rtransform!(df, :parameter = string(:parameter)) # so it is treated as categorical
    colmetadata!(df, :parameter, "label", latexstring(paramstring); style=:note) # modifies legend title

    ### misclassification rate - mean & standard deviation
    lines_mr = visual(Lines) * mapping(:T, :mr_mean; color=:parameter)
    band_mr = visual(Band; alpha=0.3) * mapping(:T, :mr_lower, :mr_upper; color=:parameter)
    lb_mr = data(df) * (band_mr + lines_mr)

    ### exact recovery - mean & confidence interval
    lines_er = visual(Lines; linestyle=:dash) * mapping(:T, :er_mean; color=:parameter)
    band_er = visual(Band; alpha=0.3) * mapping(:T, :er_lower, :er_upper; color=:parameter)
    lb_er = data(df) * (band_er + lines_er)

    ### plot
    fig, grid = draw(
        lb_mr + lb_er,
        scales(; X=(; label=L"T"), Y=(; label=""));
        axis=(; width=500, height=300, limits=(nothing, (0.0, 1.0))),
        legend=(; show=false), # remove outside legend
    )

    ### add inside legend
    ax = fig[1, 1]
    legend!(ax, grid; halign=:left, valign=:top, margin=(10, 10, 10, 10))

    return fig
end

function plot_heatmap(
    df_mean_bands::DataFrame,
    feature::Symbol,
    method::Symbol,
    clustering::Symbol;
    level::Float64=0.1,
)
    feature in [:er, :mr] || ArgumentError(
        "Only exact recovery (:er) or misclassification rate (:mr) are supported features.",
    )

    ### filter data frame for method and clustering of interest
    df = @rsubset(
        df_mean_bands, :method == string(method), :clustering == string(clustering)
    )

    ### heatmap
    hm = data(df) * visual(Heatmap) * mapping(:T, :parameter, Symbol(feature, :_mean))

    ### ranges
    tmin = minimum(unique(df.T))
    tmax = maximum(unique(df.T))
    nmax = maximum(unique(df.parameter))
    nmin = minimum(unique(df.parameter))

    ### feature specific variables
    if feature == :er
        ts = range(tmin, tmax, 100)
        ns = sqrt.(ts)
        line_lab = L"T = N^2"
        fig_title = "Exact recovery"
        color_lab = "Probability"
    else
        id_ts = @chain df begin
            @groupby(:parameter)
            @combine(:tmp = findfirst(:mr_mean .< level))
            _[!, :tmp]
        end
        ts = map(id -> isnothing(id) ? missing : unique(df.T)[id], id_ts)  # add missing values for T which are not plotted
        ns = unique(df.parameter)
        line_lab = "MMR = " * string(level)
        fig_title = "Misclassification rate"
        color_lab = "Mean"
    end

    ### line
    df_line = (T=ts, parameter=ns)
    line =
        data(df_line) * mapping(:T, :parameter) * visual(Lines; color=:red, label=line_lab)

    ### plot
    colmap = feature == :er ? :viridis : Reverse(:viridis)
    colrange = feature == :er ? (0.0, 1.0) : (0.0, 0.5)
    fig, grid = draw(
        scales(;
            X=(; label=L"T"),
            Y=(; label=L"N"),
            Color=(; label=color_lab, colorrange=colrange, colormap=colmap),
        );
        figure=(; title=fig_title, titlealign=:center),
        axis=(; width=300, height=300, limits=((tmin, tmax), (nmin, nmax))),
        legend=(; show=false), # remove outside legend
    )(
        hm + line
    )

    ### add inside legend
    ax = fig[1, 1]
    legend!(ax, grid; halign=:left, valign=:top, margin=(10, 10, 10, 10))

    return fig
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
    df = classiferrortable(
        Paramsymbol,
        Paramvec,
        default_values,
        tvec;
        multi_thread=true,
        methods=[:ag],
        clusterings=[:kmeans],
    )
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