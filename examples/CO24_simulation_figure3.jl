include("functions_classification.jl")

## Informations for reproducibility
# gitcommit : e01748726f8295d48d30155bb376f073d2f7fd34

## Default values
default_values = (N = 50,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Nsimu = Int(1e3))

# Specific values
T = Int(1e4)
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin, T, length_tvec)))

Random.seed!(1)
simulationandsave(:N, [10, 100, 200], default_values, tvec)
Random.seed!(1)
simulationandsave(:r₊, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:β, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:λ, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:p, [.1, .4, .6, .9], default_values, tvec)