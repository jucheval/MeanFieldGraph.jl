include("functions_classification.jl")
using Latexify

# Helper functions
"""
    non_optimal_couples(df, tup, metric; aside_tuples=Tuple{Symbol,Symbol}[])

Give the couples (N,T) for which the method and clustering combination given in `tup` does not have the optimal value of the metric given in `metric` (either :er or :mr), excluding the combinations given in `aside_tuples`.
"""
function non_optimal_couples(
    df::DataFrame,
    tup::Tuple{Symbol,Symbol},
    metric::Symbol;
    aside_tuples::Vector{Tuple{Symbol,Symbol}}=Tuple{Symbol,Symbol}[],
)
    df_copy = deepcopy(df) # to avoid modifying the original dataframe with @rsubset!
    f = metric == :er ? identity : (x -> -x) # to maximize ER and minimize MR

    for a_tup in aside_tuples
        @rsubset!(df_copy, !((:method, :clustering) == string.(a_tup))) # exclude the aside methods and clusterings to compare to
    end

    output = @chain df_copy begin
        transform(Symbol(metric, :_mean) => f => :target_metric) # create a column with the target metric to maximize (either ER or minus MR)
        @aside target = @rsubset(_, (:method, :clustering) == string.(tup))
        @groupby([:parameter, :T])
        @combine(:optimal_target_metric = maximum(:target_metric))
        _[target[!, :target_metric] .< _[!, :optimal_target_metric], [:parameter, :T]]
        @rename(:N = :parameter)
    end

    return output
end

round_percent(x) = round(Int, 100 * x)

## Load table
df_wide = estimatorsload("data/CO24/data_for_color_plot")
df_mean_bands = mmr_per(df_wide)

## Misclassification rate - table
df_mr = @chain df_mean_bands begin
    # value of interest is MMR ± standard error
    @rselect(
        :N = :parameter,
        :T = :T,
        :method = :method,
        :clustering = :clustering,
        :mr = string(round_percent(:mr_mean)) * " ± " * string(round_percent(:mr_std)),
    )
    # create a column with the combination of clustering and method to be able to unstack both at the same time
    @rtransform!(:col_name = string(:method, "_", :clustering))
    # unstack the table to have one column per method and clustering combination
    unstack([:N, :T], :col_name, :mr)
end

# couples are selected so that the lowest MMR (usually ag_kmeans) is around 0.02
function isselected(n, t)
    return (n, t) in [(34, 2641), (94, 7903), (142, 10534), (190, 13165), (250, 18427)]
end
latexify(@rsubset(df_mr, isselected(:N, :T)))

@chain df_mr begin
    @rsubset(isselected(:N, :T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end

## Misclassification rate - best method ?
### Most of the time, ag_threshold is amongst the lowest Mean Misclassification Rate (MMR)
### 25 couples (N, T) out of 420 where ag_threshold is not amongst the lowest MMR
size(non_optimal_couples(df_mean_bands, (:ag, :threshold), :mr))[1]
#### If we remove threshold clustering, 
#### 27 couples where ag_kmeans is not amongst the lowest MMR
size(
    non_optimal_couples(
        df_mean_bands,
        (:ag, :kmeans),
        :mr;
        aside_tuples=[(:ag, :threshold), (:sp, :threshold)],
    ),
)[1]

### When MMR is rounded up to 2 decimals, 
df_rounded = @rtransform(df_mean_bands, :mr_mean = round_percent(:mr_mean))
### only 13 couples where ag_threshold is not amongst the lowest MMR
size(non_optimal_couples(df_rounded, (:ag, :threshold), :mr))[1]
#### If we remove threshold clustering, 
#### 14 couples where ag_kmeans is not amongst the lowest MMR
size(
    non_optimal_couples(
        df_rounded, (:ag, :kmeans), :mr; aside_tuples=[(:ag, :threshold), (:sp, :threshold)]
    ),
)[1]

## Exact recovery
factor = quantile(Normal(), 0.975) / sqrt(metadata(df_mean_bands, "Number of simulations"))
df_er = @chain df_mean_bands begin
    # value of interest is PER ± standard confidence interval radius
    @rselect(
        :N = :parameter,
        :T = :T,
        :method = :method,
        :clustering = :clustering,
        :er =
            string(round_percent(:er_mean)) *
            " ± " *
            string(round_percent(factor * :er_std)),
    )
    # create a column with the combination of clustering and method to be able to unstack both at the same time
    @rtransform!(:col_name = string(:method, "_", :clustering))
    # unstack the table to have one column per method and clustering combination
    unstack([:N, :T], :col_name, :er)
end

# couples are selected so that the highest PER (usually ag_kmeans) is around 90%
function isselected(n, t)
    return (n, t) in [(34, 5272), (94, 15796), (142, 26320), (190, 36844), (250, 50000)]
end
latexify(@rsubset(df_er, isselected(:N, :T)))

@chain df_er begin
    @rsubset(isselected(:N, :T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end

## Exact recovery - best method ?
### Most of the time, ag_threshold is amongst the highest Probability of Exact Recovery (PER)
### 7 couples (N, T) out of 420 where ag_threshold is not amongst the highest PER
size(non_optimal_couples(df_mean_bands, (:ag, :threshold), :er))[1]
#### If we remove threshold clustering, 
#### 11 couples where ag_kmeans is not amongst the highest PER
size(
    non_optimal_couples(
        df_mean_bands,
        (:ag, :kmeans),
        :er;
        aside_tuples=[(:ag, :threshold), (:sp, :threshold)],
    ),
)[1]

### When PER is rounded up to 2 decimals, 
df_rounded = @rtransform(df_mean_bands, :er_mean = round_percent(:er_mean))
### only 1 couple where ag_threshold is not amongst the highest PER
size(non_optimal_couples(df_rounded, (:ag, :threshold), :er))[1]
#### If we remove threshold clustering, 
#### 2 couples where ag_kmeans is not amongst the highest PER
size(
    non_optimal_couples(
        df_rounded, (:ag, :kmeans), :er; aside_tuples=[(:ag, :threshold), (:sp, :threshold)]
    ),
)[1]