using MeanFieldGraph
using Distributions
using DataFrames
using ProgressLogging
using StatsPlots
using Profile

N = 200
r₊ = 0.6
β = 0.5
λ = 0.5
p = 0.3
μ = β * λ
model = MarkovChainModel(μ, λ, p)
m, v, w = mvw(model, r₊)

MeanFieldGraph.Φ(m, v, w, r₊)

## Typical range of parameters
# firing rate between 10 and 100 Hz -> 50 Hz
# time bin : 2 ms -> m should be close to 0.1
m = 0.1
# r₊ between 0.6 and 0.8 -> 0.7
r₊ = 0.7

N = 500
# Record time length = 20 sec -> T close to 20 sec / 2 ms = 10 000
T = Int(1e4)
Nsimu = Int(1e2)
Δ = 1
type = "absolute"
length_tvec = 100
tmin = 10
tvec = floor.(Int, collect(range(tmin, T, length_tvec)))
plow = 0.25
pup = 0.75

function βconstrained(λ, p, r₊, m)
    r₋ = 1 - r₊
    num = m * (1 - (1 - λ) * p * (r₊ - r₋)) - (1 - λ) * p * r₋
    return num / λ
end

## Compute the analytical formula of the domain (λ, p) on which βconstrained(λ, p, r₊, m) < 0
# βconstrained(λ, p, r₊, m) < 0     iff     (1-λ)p > m / (m*(r₊ - r₋) + r₋) 

for (λ, p) in Iterators.product(0.1:0.1:0.9, 0.1:0.1:0.9)
    β = βconstrained(λ, p, r₊, m)
    β < 0 && break
    println((1 - λ) * p)
    model = DiscreteHawkesModel(λ * β, λ, p)
end
