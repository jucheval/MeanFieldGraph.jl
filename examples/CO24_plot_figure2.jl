include("functions_classification.jl")

df = estimatorsload("data/CO24/data_for_color_plot.csv")
errors_and_std = proportion2errors(df)
T = errors_and_std.T
N = errors_and_std.parameter
proba = errors_and_std.proba_exact_recovery

using GLMakie
begin
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel = "T", ylabel = "N", title = "Probability of exact recovery")
    hm = Makie.heatmap!(T, N, proba, interpolate = true)
    lines!(ax, 0:100:T[end], sqrt.(0:100:T[end]), color = :red)
    Makie.xlims!(ax, tmin, T[end])
    Colorbar(fig[:, end+1], hm)
    fig
end
# line which correspond to T = 4*N^2