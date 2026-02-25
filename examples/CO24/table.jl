include("functions_classification.jl")

## Load table
df_wide = estimatorsload("data/CO24/data_for_color_plot")
df_mean_bands = mmr_per(df_wide)

### helper function
round_percent(x) = round(Int, 100 * x)

## Misclassification rate - table
df_mr = @chain df_mean_bands begin
    @transmute(
        N = parameter,
        T = T,
        method = method,
        clustering = clustering,
        mr = string(round_percent(mr_mean)) * " ± " * string(round_percent(mr_std)),
    )
    @pivot_wider(names_from = method, values_from = mr)
    @pivot_wider(names_from = clustering, values_from = ag:sp)
end

function isselected(n, t) # couples are selected so that the lowest MMR (usually ag_kmeans) is around 0.02
    return (n, t) in [(34, 2641), (94, 7903), (142, 13165), (190, 21058), (250, 31582)]
end
latexify(@filter(df_mr, isselected(N, T)))

@chain df_mr begin
    @filter(isselected(N, T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end

## Misclassification rate - best method ?
### Most of the time, the lowest Mean Misclassification Rate (MMR) is achieved for ag_kmeans
### 25 couples (N, T) out of 420 where the lowest MMR is not ag_kmeans
@chain df_mean_bands begin
    @group_by([:parameter, :T])
    @summarize(id_lowest_mmr = argmin(mr_mean))
    @ungroup
    @filter(id_lowest_mmr != 1)
end

### When MMR is rounded up to 2 decimals, only 11 couples (N, T) out of 420
@chain df_mean_bands begin
    @group_by([:parameter, :T])
    @summarize(id_lowest_mmr = argmin(round(100 * mr_mean)))
    @ungroup
    @filter(id_lowest_mmr != 1)
end

## Exact recovery
df_er = @chain df_mean_bands begin
    @transmute(
        N = parameter,
        T = T,
        method = method,
        clustering = clustering,
        er = string(round_percent(er_mean)) * " ± " * string(round_percent(er_std)),
    )
    @pivot_wider(names_from = method, values_from = er)
    @pivot_wider(names_from = clustering, values_from = ag:sp)
end

function isselected(n, t) # couples are selected so that the highest PER (usually ag_kmeans) is around 90%
    return (n, t) in [(34, 5272), (94, 15796), (142, 26320), (190, 36844), (250, 50000)]
end
latexify(@filter(df_er, isselected(N, T)))

@chain df_er begin
    @filter(isselected(N, T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end

## Exact recovery - best method ?
### Most of the time, the highest Probability of Exact Recovery (PER) is achieved for ag_kmeans
### 11 couples (N, T) out of 420 where the highest PER is not ag_kmeans
@chain df_mean_bands begin
    @group_by([:parameter, :T])
    @summarize(id_highest_mmr = argmax(er_mean))
    @ungroup
    @filter(id_highest_mmr != 1)
end

### When PER is rounded up to 2 decimals, ag_kmeans is always among the best
@chain df_mean_bands begin
    @group_by([:parameter, :T])
    @summarize(id_highest_er = argmax(round(100 * er_mean)))
    @ungroup
    @filter(id_highest_er != 1)
end