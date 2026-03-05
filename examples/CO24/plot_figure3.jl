include("functions_classification.jl")

## N varies
paramstring = "N"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
df = @rename(df, :miscl_rate_ag_kmeans = :prop_errors) # comply with new version of column names
df_mean_bands = mmr_per(df)
fig = plotclassification(df_mean_bands)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## r₊ varies
paramstring = "r₊"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
df = @rename(df, :miscl_rate_ag_kmeans = :prop_errors) # comply with new version of column names
df_mean_bands = mmr_per(df)
fig = plotclassification(df_mean_bands)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## β varies
paramstring = "β"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
df = @rename(df, :miscl_rate_ag_kmeans = :prop_errors) # comply with new version of column names
df_mean_bands = mmr_per(df)
fig = plotclassification(df_mean_bands)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## λ varies
paramstring = "λ"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
df = @rename(df, :miscl_rate_ag_kmeans = :prop_errors) # comply with new version of column names
df_mean_bands = mmr_per(df)
fig = plotclassification(df_mean_bands)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## p varies
paramstring = "p"
df = estimatorsload("data/CO24/classification_vary_" * paramstring)
df = @rename(df, :miscl_rate_ag_kmeans = :prop_errors) # comply with new version of column names
df_mean_bands = mmr_per(df)
fig = plotclassification(df_mean_bands)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)