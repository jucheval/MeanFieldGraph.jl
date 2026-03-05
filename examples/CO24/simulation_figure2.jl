include("functions_classification.jl")

### Informations for reproducibility
# gitcommit : XXX

#region Left plot
## Simulation (multi-threaded)
default_values = (N=50, r₊=0.5, β=0.5, λ=0.5, p=0.5, Nsimu=Int(1e3))

# Specific values
Nvec = 10:12:250
T = Int(5e4)
length_tvec = 20
tmin = 10
tvec = floor.(Int, collect(range(tmin, T, length_tvec)))

Random.seed!(1)
for N in Nvec ## Since it takes a long time, a CSV is saved for each N
    df = classiferrortable(:N, [N], default_values, tvec; multi_thread=true)

    # Save
    CSV.write("data/CO24/data_for_color_plot_N=" * string(N) * ".csv", df)
    open("data/CO24/data_for_color_plot_N=" * string(N) * ".toml", "w") do io
        print(io, meta2toml(df))
    end
end

## Load and merge individual CSVs into a large one and save it
begin
    df = DataFrame()
    # Load and merge
    for N in Nvec
        df_partial = CSV.read(
            "data/CO24/data_for_color_plot_N=" * string(N) * ".csv", DataFrame
        )
        open("data/CO24/data_for_color_plot_N=" * string(N) * ".toml") do io
            toml2meta!(df_partial, io)
        end
        append!(df, df_partial)
    end
    metadatacomplete!(df, :N, default_values)

    # Save
    CSV.write("data/CO24/data_for_color_plot.csv", df)
    open("data/CO24/data_for_color_plot.toml", "w") do io
        print(io, meta2toml(df))
    end
end
#endregion

#region Right plot
## Simulation (multi-threaded)
default_values = (N=50, r₊=0.5, β=0.5, λ=0.5, p=0.5, Nsimu=Int(1e3))

# Specific values
Nvec = 10:12:250
T = Int(5e3)
length_tvec = 100
tmin = 10
tvec = floor.(Int, collect(range(tmin, T, length_tvec)))

Random.seed!(1)
for N in Nvec ## Since it takes a long time, a CSV is saved for each N
    df = classiferrortable(:N, [N], default_values, tvec; multi_thread=true, methods=[:ag]) # restrict to aggregated method to save computational time

    # Save
    CSV.write("data/CO24/data_for_color_plot_mr_N=" * string(N) * ".csv", df)
    open("data/CO24/data_for_color_plot_mr_N=" * string(N) * ".toml", "w") do io
        print(io, meta2toml(df))
    end
end

## Load and merge individual CSVs into a large one and save it
begin
    df = DataFrame()
    # Load and merge
    for N in Nvec
        df_partial = CSV.read(
            "data/CO24/data_for_color_plot_mr_N=" * string(N) * ".csv", DataFrame
        )
        open("data/CO24/data_for_color_plot_mr_N=" * string(N) * ".toml") do io
            toml2meta!(df_partial, io)
        end
        append!(df, df_partial)
    end
    metadatacomplete!(df, :N, default_values)

    # Save
    CSV.write("data/CO24/data_for_color_plot_mr.csv", df)
    open("data/CO24/data_for_color_plot_mr.toml", "w") do io
        print(io, meta2toml(df))
    end
end
#endregion