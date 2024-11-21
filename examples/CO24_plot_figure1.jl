include("functions_classification.jl")

# N varies
paramstring = "N"
df = estimatorsload("data/CO24/classification_vary_"*paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)

# r₊ varies
paramstring = "r₊"
df = estimatorsload("data/CO24/classification_vary_"*paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)

# β varies
paramstring = "β"
df = estimatorsload("data/CO24/classification_vary_"*paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)

# λ varies
paramstring = "λ"
df = estimatorsload("data/CO24/classification_vary_"*paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)

# p varies
paramstring = "p"
df = estimatorsload("data/CO24/classification_vary_"*paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)