"""
    classification(data::DiscreteTimeData)

Estimates the two underlying communities (one excitatory and one inhibitory) from the data set `data`. It returns a `Vector{Bool}` where the `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.

# Keyword arguments
    - `clustering::Symbol`: the clustering method applied to the estimated covariance vector. Valid choices are: `:kmeans` (the default) and `:hclust`.
"""
function classification(data::DiscreteTimeData; clustering::Symbol=:kmeans)::Vector{Bool}
    σ̂ = covariance_vector(data)

    if clustering == :kmeans
        initialisation = [argmin(σ̂), argmax(σ̂)]
        output = cluster2bool(kmeans(transpose(σ̂), 2; init=initialisation))
    elseif clustering == :hclust
        distances = [abs(σ̂[i] - σ̂[j]) for i in eachindex(σ̂), j in eachindex(σ̂)]
        ct = cutree(hclust(distances); k=2)
        id_excitatory = ct[argmax(σ̂)]
        output = ct .== id_excitatory
    else
        throw(ArgumentError("Unsupported clustering method $clustering"))
    end

    return output
end

function covariance_vector(data::DiscreteTimeData)::Vector{Float64}
    X = data.X
    N, T = size(data)
    Z = sum(X; dims=2)
    ΣX = sum(X; dims=1)

    s = @views(X[:, 1:(end - 1)] * ΣX[2:end])
    output = s / (T - 1) - sum(Z) * Z / T^2
    # It is not written exactly like the definition in the paper but it is an equivalent formula.

    return output[:, 1]
end

# Auxiliary functions

function cluster2bool(R::ClusteringResult)::Vector{Bool}
    output = Vector{Bool}(undef, sum(counts(R)))
    check = R.centers[1] < R.centers[2]

    output[assignments(R) .== 1] .= !check
    output[assignments(R) .== 2] .= check

    return output
end
