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

begin # Plot
    Plots.heatmap(unique(T), unique(N), transpose(reshape(proba, (20,21))))
    plot!(0:100:T[end], sqrt.(0:100:T[end]), color = :green1, label=L"T = N^2")
    xlims!(T[1], T[end])
    ylims!(N[1], N[end])
    xlabel!(L"T")
    ylabel!(L"N")
end
savefig("plots/CO24/colorplot.pdf")