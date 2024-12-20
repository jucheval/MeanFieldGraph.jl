include("functions_estimation.jl")

# Values for the quartiles for the enveloppe
plow = 0.25
pup = 0.75
# How to add ribbons to the plots ? 
# Replace 
# plotestimatorserrors(df, df_inf)
# by
# plotestimatorserrors(df, df_inf; quantiles=(plow, pup))
# below

# N varies
paramstring = "N"
df, df_inf = estimatorsload("data/CLO24/estimators_vary_" * paramstring * "_delta_1")
plotestimatorserrors(df, df_inf)

# r₊ varies
paramstring = "r₊"
df, df_inf = estimatorsload("data/CLO24/estimators_vary_" * paramstring * "_delta_1")
plotestimatorserrors(df, df_inf)

# β varies
paramstring = "β"
df, df_inf = estimatorsload("data/CLO24/estimators_vary_" * paramstring * "_delta_1")
plotestimatorserrors(df, df_inf)

# λ varies
paramstring = "λ"
df, df_inf = estimatorsload("data/CLO24/estimators_vary_" * paramstring * "_delta_1")
plotestimatorserrors(df, df_inf)

# p varies
paramstring = "p"
df, df_inf = estimatorsload("data/CLO24/estimators_vary_" * paramstring * "_delta_1")
plotestimatorserrors(df, df_inf)
