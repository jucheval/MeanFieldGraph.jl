include("functions_classification.jl")

"""
    plot_vary(s::String)

Helper function to plot the performance of the `ag_kmeans` method as a function of a varying parameter `s`.
"""
function plot_vary(s::String)
    df = estimatorsload("data/CO24/classification_vary_" * s)
    df_mean_bands = mmr_per(df)
    @rsubset!(df_mean_bands, :method .== "ag" && :clustering .== "kmeans")  # select ag_kmeans only
    return plotclassification(df_mean_bands)
end

## N varies
paramstring = "N"
fig = plot_vary(paramstring)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## r₊ varies
paramstring = "r₊"
fig = plot_vary(paramstring)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## β varies
paramstring = "β"
fig = plot_vary(paramstring)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## λ varies
paramstring = "λ"
fig = plot_vary(paramstring)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)

## p varies
paramstring = "p"
fig = plot_vary(paramstring)
save("plots/CO24/classification_vary_" * paramstring * ".pdf", fig)