"""
    classification(data::DiscreteTimeData)

Estimates the two underlying communities (one excitatory and one inhibitory) from the data set `data`. It returns a `Vector{Bool}` where the `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
"""
function classification(data::DiscreteTimeData)::Vector{Bool}
    N, T = size(data)

    σ̂ = covariance_vector(data)
    τ = sortperm(σ̂)
    D = diff(σ̂[τ])
    k̂ = argmax(D)

    output = ones(Bool, N)
    output[τ[1:k̂]] .= false

    return output
end

function covariance_vector(data::DiscreteTimeData)::Vector{Float64}
    X = data.X
    N, T = size(data)
    Z = sum(X; dims=2)
    ΣX = sum(X; dims=1)

    X = X[:,1:end-1]
    ΣX = ΣX[2:end]
    s = X * ΣX
    output = s/(T-1) - sum(Z)*Z/T^2
    # It is not written exactly like the definition in the paper but it is an equivalent formula.
    
    return output[:,1]
end