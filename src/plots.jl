import Plots.plot

"""
    plot(data::DiscreteTimeData)

Plot a black and white image corresponding to `data`. Black pixels correpond to the value `1`. The x-axis correspond to time and y-axis correpond to components.
"""
function plot(data::DiscreteTimeData)
    N, T = size(data)
    heatmap(data.X, color=palette(:binary), xticks = floor.(Int,range(1,T,10)), yticks = floor.(Int,range(1,N,10)), colorbar=:false, tick_direction=:none)
end