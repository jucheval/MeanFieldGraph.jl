"""
    classification(data::DiscreteTimeData)

Estimates the two underlying communities (one excitatory and one inhibitory) from the data set `data`. It returns a `Vector{Bool}` where the `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.

# Keyword arguments
    - `method::Symbol`: the method applied to estimate the covariance vector σ. Valid choices are: `:aggregated` (the default) and `:spectral`.
    - `clustering::Symbol`: the clustering method applied to the estimated covariance vector. Valid choices are: `:kmeans` (the default) and `:hclust`.
"""
function classification(
    data::DiscreteTimeData; method::Symbol=:aggregated, clustering::Symbol=:kmeans
)::Vector{Bool}

    # Estimation of the covariance vector σ
    if method == :aggregated
        σ̂ = covariance_vector(data)
    elseif method == :spectral
        Σ̂ = covariance_matrix(data)
        _, vecs = eigsolve(transpose(Σ̂) * Σ̂) # faster than full SVD
        σ̂ = vecs[1]
        # FIXME : how to chose the sign of σ̂sp ?
    else
        throw(ArgumentError("Unsupported method $method"))
    end

    # Clustering based on the estimated σ̂
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

function covariance_matrix(data::DiscreteTimeData)::Matrix{Float64}
    X = data.X
    N, T = size(data)
    Z = sum(X; dims=2)

    s = zeros((N, N))
    for t in 1:(T - 1)
        s += @views(X[:, t + 1] * transpose(X[:, t]))
    end
    output = s / (T - 1) - Z * transpose(Z) / T^2

    return output
end

# Auxiliary functions

function cluster2bool(R::ClusteringResult)::Vector{Bool}
    output = Vector{Bool}(undef, sum(counts(R)))
    check = R.centers[1] < R.centers[2]

    output[assignments(R) .== 1] .= !check
    output[assignments(R) .== 2] .= check

    return output
end
