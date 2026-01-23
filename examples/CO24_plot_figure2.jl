include("functions_classification.jl")

## Left plot
df_wide = estimatorsload("data/CO24/data_for_color_plot")
df_mean_bands = mmr_per(df_wide)

fig = plot_heatmap(df_mean_bands, :er, :ag, :kmeans)

save("plots/CO24/heatmap_er.pdf", fig)

## Right plot
df_wide = estimatorsload("data/CO24/data_for_color_plot_mr")
df_mean_bands = mmr_per(df_wide)

fig = plot_heatmap(df_mean_bands, :mr, :ag, :kmeans)

save("plots/CO24/heatmap_mr.png", fig)