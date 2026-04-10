include("functions_classification.jl")

## Informations for reproducibility
# gitcommit : d9cdda27e086151d37b742fa72e4c182a6ebfe5c

## Default values
default_values = (N=50, r₊=0.5, β=0.5, λ=0.5, p=0.5, Nsimu=Int(1e3))

# Specific values
T = Int(1e4)
length_tvec = 100
tmin = 10
tvec = floor.(Int, collect(range(tmin, T, length_tvec)))

# method = :ag, all clustering
Random.seed!(1)
simulationandsave(:N, [10, 100, 200], default_values, tvec)
Random.seed!(1)
simulationandsave(:r₊, [0.1, 0.4, 0.6, 0.9], default_values, tvec)
Random.seed!(1)
simulationandsave(:β, [0.1, 0.4, 0.6, 0.9], default_values, tvec)
Random.seed!(1)
simulationandsave(:λ, [0.1, 0.4, 0.6, 0.9], default_values, tvec)
Random.seed!(1)
simulationandsave(:p, [0.1, 0.4, 0.6, 0.9], default_values, tvec)
