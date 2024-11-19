include("functions_estimation.jl")
Random.seed!(1)

## Informations for reproducibility

## Default values
N = 500
r₊ = .5
β = .5
λ = .5
p = .5
Δ = 1
Nsimu = Int(1e3)

# Specific values
model = MarkovChainModel(β*λ,λ,p)
Nsimu = 10
T = Int(1e5)
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))
Δvec = [0,1]

# Simulation
df, df_inf = estimatorstable(model, N, r₊, Nsimu, tvec, Δvec)

# Plots
df_complete_inf = extractw_computeμλp(df_inf, 1)
estimators2errors!(df_complete_inf, model, r₊)
# Left plot
df_left = extractw_computeμλp(df, 1)
estimators2errors!(df_left, model, r₊)
df_median = columnquantile(df_left, .5)
df_median_inf = columnquantile(df_complete_inf, .5)

plot(yaxis=:log, legend_columns=3)
for id in 1:6
    plot!(df_median.T, df_median[:,id+1], color=id, label=names(df_median)[id+1])
    scatter!(T*ones(size(df_median_inf)[1]), df_median_inf[:,id], color=id, marker = :o, label=false)
end
xlabel!("T")
ylims!(1e-4, 1)
ylabel!("absolute error")
title!("Estimation error for all parameters, Δ="*string(Δ))