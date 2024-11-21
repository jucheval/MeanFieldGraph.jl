include("functions_estimation.jl")
using Logging
Random.seed!(1)

## Informations for reproducibility

## Default values
default_values = (N = 500,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Δ = 1,
Nsimu = Int(1e3))

# Specific values
T = Int(1e3)
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))

disable_logging(LogLevel(-1001))  # Enables debug info

# Simulation and saving
simulationandsave(:N, [100, 1000, 2000], default_values, tvec)
simulationandsave(:r₊, [.1, .4, .6, .9], default_values, tvec)
simulationandsave(:β, [.1, .4, .6, .9], default_values, tvec)
simulationandsave(:λ, [.1, .4, .6, .9], default_values, tvec)
simulationandsave(:p, [.1, .4, .6, .9], default_values, tvec)