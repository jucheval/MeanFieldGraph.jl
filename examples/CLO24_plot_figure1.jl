include("functions_estimation.jl")
begin # Load 
    df = CSV.read("data/CLO24/estimators_vary_Δ.csv", DataFrame)
    open("data/CLO24/estimators_vary_Δ.toml") do io
        toml2meta!(df, io)
    end
    df_inf = CSV.read("data/CLO24/estimators_vary_Δ_inf.csv", DataFrame)
    open("data/CLO24/estimators_vary_Δ_inf.toml") do io
        toml2meta!(df_inf, io)
    end
end;

begin # Compute absolute errors and medians
    errors, errors_inf = estimators2errors.((df, df_inf))
    q50, q50_inf = columnquantile.((errors, errors_inf), 0.5)
end;

begin # Left plot 
    selection = (q50.parameter .== 0)
    plot(; yaxis=:log, legend_columns=3)
    for id in 1:6
        plot!(q50[selection, 1], q50[selection, id + 2]; color=id, label=names(q50)[id + 2])
        scatter!((q50[end, 1], q50_inf[1, id + 1]); color=id, marker=:o, label=false)
    end
    xlabel!("T")
    ylims!(1e-4, 1)
    ylabel!("absolute error")
    title!("Estimation error for all parameters, Δ=NaN")
end;

begin # Right plot 
    selection = (q50.parameter .== 1)
    plot(; yaxis=:log, legend_columns=3)
    for id in 1:6
        plot!(q50[selection, 1], q50[selection, id + 2]; color=id, label=names(q50)[id + 2])
        scatter!((q50[end, 1], q50_inf[1, id + 1]); color=id, marker=:o, label=false)
    end
    xlabel!("T")
    ylims!(1e-4, 1)
    ylabel!("absolute error")
    title!("Estimation error for all parameters, Δ=1")
end;
