"""
    classification(data::DiscreteTimeData)

Estimates the two underlying communities (one excitatory and one inhibitory) from the data set `data`. It returns a `Vector{Bool}` where the `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.

# Keyword arguments
    - `method::Symbol`: the method applied to estimate the covariance vector σ. Valid choices are: `:aggregated` (the default) and `:spectral`.
    - `clustering::Symbol`: the clustering method applied to the estimated covariance vector. Valid choices are: `:kmeans` (the default), `:threshold`, and the *linkage* choices for the `hclust` function (`:single`, `:average`, `:complete`, `:ward`).
"""
function classification(
    data::DiscreteTimeData; method::Symbol=:aggregated, clustering::Symbol=:kmeans
)::Vector{Bool}

    # Estimation of the covariance vector σ
    if method == :aggregated
        σ̂ = covariance_vector(data)
    elseif method == :spectral
        # compute the leading singular vector of the covariance matrix
        Σ̂ = covariance_matrix(data)
        _, vecs = eigsolve(transpose(Σ̂) * Σ̂) # faster than full SVD
        v̌ = vecs[1]

        # sign disambiguation
        σ̂_ag = sum(Σ̂; dims=1)[1, :]
        m_v̌ = mean(v̌)
        P̌ = v̌ .>= m_v̌
        σ̌₊ = sum(σ̂_ag[P̌]) / sum(P̌)
        σ̌₋ = sum(σ̂_ag[.!P̌]) / sum(.!P̌)

        σ̂ = σ̌₊ >= σ̌₋ ? v̌ : -v̌
    else
        throw(ArgumentError("Unsupported method $method"))
    end

    # Clustering based on the estimator σ̂
    if clustering == :kmeans
        initialisation = [argmin(σ̂), argmax(σ̂)]
        output = cluster2bool(kmeans(transpose(σ̂), 2; init=initialisation))
    elseif clustering == :threshold
        threshold = mean(σ̂)
        output = σ̂ .>= threshold
    elseif clustering in (:single, :average, :complete, :ward)
        distances = pairwise_abs_distances(σ̂)
        ct = cutree(hclust(distances; linkage=clustering); k=2)
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

    s = @views(X[:, 2:end] * transpose(X[:, 1:(end - 1)]))
    output = s / (T - 1) - Z * transpose(Z) / T^2

    return output
end

# Auxiliary functions

function pairwise_abs_distances(v::Vector{<:Real})::Matrix{Float64}
    n = length(v)
    output = Matrix{Float64}(undef, n, n)
    @inbounds for i in 1:n
        output[i, i] = 0.0
        @inbounds for j in (i + 1):n
            d = abs(v[i] - v[j])
            output[i, j] = d
            output[j, i] = d
        end
    end
    return output
end

function cluster2bool(R::ClusteringResult)::Vector{Bool}
    output = Vector{Bool}(undef, sum(counts(R)))
    check = R.centers[1] < R.centers[2]

    output[assignments(R) .== 1] .= !check
    output[assignments(R) .== 2] .= check

    return output
end
