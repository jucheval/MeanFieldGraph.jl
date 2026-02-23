include("functions_classification.jl")

## Load table
df_wide = estimatorsload("data/CO24/data_for_color_plot")
df_mean_bands = mmr_per(df_wide)

### helper function
round_percent(x) = round(Int, 100 * x)

## Misclassification rate
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

function isselected(n, t) # couples are selected so that the lowest MMR (usually kmeans_ag) is around 0.02
    return (n, t) in [(34, 2641), (94, 7903), (142, 13165), (190, 21058), (250, 31582)]
end
latexify(@filter(df_mr, isselected(N, T)))

@chain df_mr begin
    @filter(isselected(N, T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
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

function isselected(n, t) # couples are selected so that the highest PER (usually kmeans_ag) is around 90%
    return (n, t) in [(34, 5272), (94, 15796), (142, 26320), (190, 36844), (250, 50000)]
end
latexify(@filter(df_er, isselected(N, T)))

@chain df_er begin
    @filter(isselected(N, T))
    latexify(; env=:table, booktabs=true, latex=true)
    print()
end