include("functions_classification.jl")

## Informations for reproducibility
# gitcommit : e01748726f8295d48d30155bb376f073d2f7fd34

## Default values
default_values = (N=250, r₊=0.5, β=0.5, λ=0.5, p=0.5, Nsimu=Int(1e3))

# Specific values
Nvec = 10:12:250
T = Int(5e4)
length_tvec = 20
tmin = 10
tvec = floor.(Int, collect(range(tmin, T, length_tvec)))

Random.seed!(1)
begin # Simulation
    df = classiferrortable(:N, Nvec, default_values, tvec)

    # Save
    CSV.write("data/CO24/data_for_color_plot.csv", df)
    open("data/CO24/data_for_color_plot.toml", "w") do io
        print(io, meta2toml(df))
    end
end;
