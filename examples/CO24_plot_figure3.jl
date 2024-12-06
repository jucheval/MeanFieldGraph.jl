include("functions_classification.jl")

# N varies
paramstring = "N"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)
savefig("plots/CO24/classification_vary_" * paramstring * ".pdf")

# r₊ varies
paramstring = "r₊"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)
savefig("plots/CO24/classification_vary_" * paramstring * ".pdf")

# β varies
paramstring = "β"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)
savefig("plots/CO24/classification_vary_" * paramstring * ".pdf")

# λ varies
paramstring = "λ"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)
savefig("plots/CO24/classification_vary_" * paramstring * ".pdf")

# p varies
paramstring = "p"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
errors_and_std = proportion2errors(df)
plotclassification(errors_and_std)
savefig("plots/CO24/classification_vary_" * paramstring * ".pdf")
