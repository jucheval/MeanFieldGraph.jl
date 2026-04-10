include("functions_classification.jl")

"""
    plot_vary(s::String)

Helper function to plot the performance of the aggregated method with clustering procedure `clustering` as a function of a varying parameter `s`.
"""
function plot_vary(s::String, clustering::Symbol=:kmeans)
    df = estimatorsload("data/CO24/classification_vary_" * s)
    df_mean_bands = mmr_per(df)
    @rsubset!(df_mean_bands, :method .== "ag" && :clustering .== string(clustering))  # select aggregated with the specified clustering only
    return plotclassification(df_mean_bands)
end

## Figure 2: performance of the aggregated method with k-means clustering as a function of the varying parameter
### N varies
paramstring = "N"
fig = plot_vary(paramstring, :kmeans)
save("plots/CO24/classification_vary_" * paramstring * "_kmeans.pdf", fig)
### r₊ varies
paramstring = "r₊"
fig = plot_vary(paramstring, :kmeans)
save("plots/CO24/classification_vary_" * paramstring * "_kmeans.pdf", fig)
### β varies
paramstring = "β"
fig = plot_vary(paramstring, :kmeans)
save("plots/CO24/classification_vary_" * paramstring * "_kmeans.pdf", fig)
### λ varies
paramstring = "λ"
fig = plot_vary(paramstring, :kmeans)
save("plots/CO24/classification_vary_" * paramstring * "_kmeans.pdf", fig)
### p varies
paramstring = "p"
fig = plot_vary(paramstring, :kmeans)
save("plots/CO24/classification_vary_" * paramstring * "_kmeans.pdf", fig)

## Figure 3: performance of the aggregated method with mean threshold clustering as a function of the varying parameter
### N varies
paramstring = "N"
fig = plot_vary(paramstring, :threshold)
save("plots/CO24/classification_vary_" * paramstring * "_threshold.pdf", fig)
### r₊ varies
paramstring = "r₊"
fig = plot_vary(paramstring, :threshold)
save("plots/CO24/classification_vary_" * paramstring * "_threshold.pdf", fig)
### β varies
paramstring = "β"
fig = plot_vary(paramstring, :threshold)
save("plots/CO24/classification_vary_" * paramstring * "_threshold.pdf", fig)
### λ varies
paramstring = "λ"
fig = plot_vary(paramstring, :threshold)
save("plots/CO24/classification_vary_" * paramstring * "_threshold.pdf", fig)
### p varies
paramstring = "p"
fig = plot_vary(paramstring, :threshold)
save("plots/CO24/classification_vary_" * paramstring * "_threshold.pdf", fig)