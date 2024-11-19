include("functions_classification.jl")
Random.seed!(1)

## Informations for reproducibility

## Default values
default_values = (N = 250,
r₊ = .5,
β = .5,
λ = .5,
p = .5,
Nsimu = 100)

# Specific values
T = Int(1e4)
length_tvec = 10
tmin = 1000
tvec = floor.(Int,collect(range(tmin, T, length_tvec)))

# N varies
paramvec = [100, 200, 500]
errors_and_std = proportion2errors(classiferrortable(:N, paramvec, default_values, tvec))
plotclassification(errors_and_std)

# r₊ varies
paramvec = [.1, .4, .6, .9]
errors_and_std = proportion2errors(classiferrortable(:r₊, paramvec, default_values, tvec))
plotclassification(errors_and_std)

# λ varies
paramvec = [.1, .4, .6, .9]
errors_and_std = proportion2errors(classiferrortable(:λ, paramvec, default_values, tvec))
plotclassification(errors_and_std)

# p varies
paramvec = [.1, .4, .6, .9]
errors_and_std = proportion2errors(classiferrortable(:p, paramvec, default_values, tvec))
plotclassification(errors_and_std)

# β varies
paramvec = [.1, .4, .6, .9]
errors_and_std = proportion2errors(classiferrortable(:β, paramvec, default_values, tvec))
plotclassification(errors_and_std)