"""
    plot(data::DiscreteTimeData)

Plot a black and white image corresponding to `data`. Black pixels correspond to
the value `1`. The x-axis corresponds to time and y-axis corresponds to
components.

The concrete plotting implementation is provided by a package extension when
`Plots` is available.
"""
function plot end
