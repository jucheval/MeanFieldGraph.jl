"""
    classification(data::DiscreteTimeData)

Estimates the two underlying communities (one excitatory and one inhibitory) from the data set `data`. It returns a `Vector{Bool}` where the `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
"""
function classification(data::DiscreteTimeData)::Vector{Bool}
    N, T = size(data)
    σ̂ = covariance_vector(data)
    initialisation = [argmin(σ̂), argmax(σ̂)]
    output = cluster2bool(kmeans(transpose(σ̂), 2, init=initialisation))
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

# Auxiliary functions

function cluster2bool(R::KmeansResult)::Vector{Bool}
    output = Vector{Bool}(undef, sum(counts(R)))
    check = R.centers[1] < R.centers[2]

    output[assignments(R) .== 1] .= !check
    output[assignments(R) .== 2] .= check

    return output
end