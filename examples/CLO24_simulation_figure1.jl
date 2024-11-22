include("functions_estimation.jl")
Random.seed!(1)

## Informations for reproducibility
# This file was created after the publication of https://hal.science/hal-04609972v1.
# Hence, the data created via this file and the associated plots are not exactly those appearing in the article.
# Nevertheless, they are rather similar.

## Default values
default_values = (N = 500,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Δ = 1,
Nsimu = Int(1e2))

# Specific values
T = Int(1e5)
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))


begin   # Simulation 
    Δvec = [0,1]
    df, df_inf = estimatorstable(:Δ, Δvec, default_values, tvec)

    paramstring = metadata(df, "Varying parameter")
    begin   # Save
        CSV.write("data/estimators_vary_Δ.csv", df)
        open("data/estimators_vary_Δ.toml", "w") do io
            print(io, meta2toml(df))
        end
    
        CSV.write("data/estimators_vary_Δ_inf.csv", df_inf)
        open("data/estimators_vary_Δ_inf.toml", "w") do io
            print(io, meta2toml(df_inf))
        end
    end;
end;