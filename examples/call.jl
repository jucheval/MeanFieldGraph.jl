using MeanFieldGraph
using Distributions
using DataFrames
using ProgressLogging
using StatsPlots
using Profile


N = 200
r₊ = .6
β = .5
λ = .5
p = .3
μ = β*λ
model = MarkovChainModel(μ,λ,p)
m,v,w = mvw(model, r₊)

MeanFieldGraph.Φ(m,v,w,r₊)


# The naive classification of the neurons by the difference between the activity after a 0 and after a 1,
# can work
excitatory = MeanFieldGraph.N2excitatory(N, r₊)
T = Int(1e6)
data = rand(model, excitatory, T)
begin
    id = 145
    time0 = pushfirst!(data.X[id, 1:(end-1)] .== 0, false)
    mean0 = mean(data.X[:, time0])
    time1 = pushfirst!(data.X[id, 1:(end-1)] .== 1, false)
    mean1 = mean(data.X[:, time1])
    println(mean0 - mean1)
end


## Typical range of parameters
# firing rate between 10 and 100 Hz -> 50 Hz
# time bin : 2 ms -> m should be close to 0.1
m=.1
# r₊ between 0.6 and 0.8 -> 0.7
r₊ = .7

N = 500
# Record time length = 20 sec -> T close to 20 sec / 2 ms = 10 000
T = Int(1e4)
Nsimu = Int(1e2)
Δ = 1
type = "absolute"
length_tvec = 100
tmin = 10
tvec = floor.(Int,collect(range(tmin,T,length_tvec)))
plow = .25
pup = .75

function βconstrained(λ, p, r₊, m)
    r₋ = 1-r₊
    num = m*(1 - (1-λ)*p*(r₊ - r₋)) - (1-λ)*p*r₋
    return num/λ
end

## Compute the analytical formula of the domain (λ, p) on which βconstrained(λ, p, r₊, m) < 0
# βconstrained(λ, p, r₊, m) < 0     iff     (1-λ)p > m / (m*(r₊ - r₋) + r₋) 

for (λ, p) in Iterators.product(.1:.1:.9,.1:.1:.9)
    β = βconstrained(λ, p, r₊, m)
    β < 0 && break
    println((1-λ)*p)
    model = DiscreteHawkesModel(λ*β, λ, p)
    # Refaire la fonction errortables pour y mettre uniquement les estimateurs m,v,w pour chaque simulation.
    # Tout écrire avec les seed de Random plutôt que de faire des CSV
    #df, df_inf = errortables(:β, type, [0], Nsimu, tvec, Δ)
end