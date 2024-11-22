include("functions_estimation.jl")
using Logging


## Informations for reproducibility
# This file was created after the publication of https://hal.science/hal-04609972v1.
# Hence, the data created via this file and the associated plots are not exactly those appearing in the article.
# Nevertheless, they are rather similar.

# The time consuming part are the five calls of the function "simulationandsave".
# The random seed is reset after each call so that they can be applied separately.

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
Random.seed!(1)
simulationandsave(:N, [100, 1000, 2000], default_values, tvec)
Random.seed!(1)
simulationandsave(:r₊, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:β, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:λ, [.1, .4, .6, .9], default_values, tvec)
Random.seed!(1)
simulationandsave(:p, [.1, .4, .6, .9], default_values, tvec)