include("functions_classification.jl")
Random.seed!(1)

## Default values
default_values = (N = 250,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Nsimu = 50)

# Specific values
T = Int(1e5)
length_tvec = 10
tmin = 1000
tvec = floor.(Int,collect(range(tmin, T, length_tvec)))

# N varies
paramvec = [100, 200, 500]
df = classiferrortables(:N, paramvec, default_values, tvec)
plotclassification(df)

# r₊ varies
paramvec = [.1, .4, .6, .9]
df = classiferrortables(:r₊, paramvec, default_values, tvec)
plotclassification(df)

# λ varies
paramvec = [.1, .4, .6, .9]
df = classiferrortables(:λ, paramvec, default_values, tvec)
plotclassification(df)

# p varies
paramvec = [.1, .4, .6, .9]
df = classiferrortables(:p, paramvec, default_values, tvec)
plotclassification(df)

# β varies
paramvec = [.1, .4, .6, .9]
df = classiferrortables(:β, paramvec, default_values, tvec)
plotclassification(df)