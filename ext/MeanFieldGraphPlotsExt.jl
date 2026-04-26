module MeanFieldGraphPlotsExt

using MeanFieldGraph
using Plots: heatmap, palette
import MeanFieldGraph: plot

function plot(data::MeanFieldGraph.DiscreteTimeData)
    N, T = size(data)
    return heatmap(
        data.X;
        color=palette(:binary),
        xticks=floor.(Int, range(1, T, 10)),
        yticks=floor.(Int, range(1, N, 10)),
        colorbar=:false,
        tick_direction=:none,
    )
end

end
