## Adapt the content of plots_p.jl with the mindset of classification_curves.jl
include("functions_estimation.jl")
using Logging
using CSV
using TableMetadataTools
Random.seed!(1)

## Informations for reproducibility

## Default values
default_values = (N = 500,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Δ = 1,
Nsimu = Int(1e1))

# Specific values
T = Int(1e3)
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))
plow = .25
pup = .75

type = "absolute"

disable_logging(LogLevel(-1001))  # Enables debug info

begin   # Simulation 
    # Change the values and the symbol of the parameter
    paramvec = [100, 1000, 2000]
    df, df_inf = estimatorstable(:N, paramvec, default_values, tvec)

    paramstring = metadata(df, "Varying parameter")
    begin   # Save
        CSV.write("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*".csv", df)
        open("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*".toml", "w") do io
            print(io, meta2toml(df))
        end
    
        CSV.write("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.csv", df_inf)
        open("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.toml", "w") do io
            print(io, meta2toml(df_inf))
        end
    end;
end;

# β varies
paramvec = [.1, .6]
df, df_inf = estimatorstable(:β, paramvec, default_values, tvec)
errors, errors_inf = estimators2errors.((df, df_inf))
q50, q50_inf = columnquantile.((errors, errors_inf), .5)
qlow, qlow_inf = columnquantile.((errors, errors_inf), plow)
qup, qup_inf = columnquantile.((errors, errors_inf), pup)
plotestimators((qlow, q50, qup), (qlow_inf, q50_inf, qup_inf))


begin # Load 
    paramstring = "N"
    df = CSV.read("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*".csv", DataFrame)
    open("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*".toml") do io
        toml2meta!(df, io)
    end
    df_inf = CSV.read("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.csv", DataFrame)
    open("../discrete-hawkes/code/data/estimators_vary_"*paramstring*"_delta_"*string(Δ)*"_inf.toml") do io
        toml2meta!(df_inf, io)
    end
end