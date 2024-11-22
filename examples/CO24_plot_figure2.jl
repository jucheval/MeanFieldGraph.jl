include("functions_classification.jl")

begin # Load
    df = CSV.read("data/CO24/data_for_color_plot.csv", DataFrame)
    open("data/CO24/data_for_color_plot.toml") do io
        toml2meta!(df, io)
    end
end;
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
    Makie.xlims!(ax, T[1], T[end])
    Makie.ylims!(ax, N[1], N[end])
    Colorbar(fig[:, end+1], hm)
    fig
end
# line which correspond to T = 4*N^2