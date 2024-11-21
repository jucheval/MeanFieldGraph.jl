include("functions_classification.jl")
Random.seed!(1)

## Informations for reproducibility

## Default values
default_values = (N = 250,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Nsimu = Int(1e3))

# Specific values
Nvec = 10:12:250
T = Int(5e4)
length_tvec = 20
tmin = 10
tvec = floor.(Int,collect(range(tmin, T, length_tvec)))

begin # Simulation
    df = classiferrortable(:N, Nvec, default_values, tvec)

    # Save
    CSV.write("data/CO24/data_for_color_plot.csv", df)
    open("data/CO24/data_for_color_plot.toml", "w") do io
        print(io, meta2toml(df))
    end
end;