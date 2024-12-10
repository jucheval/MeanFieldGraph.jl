## Structures

"""
    MarkovChainModel(μ,λ,p)

A *Markov Chain model* with parameters `μ`, `λ` and `p`.

# Arguments
- `μ`: Spontaneous event probability parameter
- `λ`: Interaction parameter
- `p`: Edge probability

# External links
* [Chevallier, Löcherbach, Ost](https://arxiv.org/abs/2406.07066)
"""
struct MarkovChainModel
    μ::Float64
    λ::Float64
    p::Float64
    function MarkovChainModel(μ::Float64, λ::Float64, p::Float64)
        if (λ < 0) || (λ > 1)
            throw(ArgumentError("λ is $λ but must be in the range [0, 1]"))
        end
        if (μ < 0) || (μ > λ)
            throw(ArgumentError("μ is $μ but must be in the range [0, $λ]"))
        end
        if (p < 0) || (p > 1)
            throw(ArgumentError("p is $p but must be in the range [0, 1]"))
        end
        return new(μ, λ, p)
    end
end

# constructors
MarkovChainModel(μ::Real, λ::Real, p::Real) = MarkovChainModel(Float64.((μ, λ, p))...)

# methods

# print
function Base.show(io::IO, ::MIME"text/plain", model::MarkovChainModel)
    return print(
        io,
        "μ = $(model.μ), λ = $(model.λ). The θ matrix is not specified (but p = $(model.p)).",
    )
end

"""
    MarkovChainConnectivity(model,θ)

A *Markov Chain model* with parameters given by `model` and connectivity matrix given by `θ`.

# Arguments
- `model::MarkovChainModel`: contains the parameters `μ`, `λ` and `p`.
- `θ::Matrix{Bool}`: `θ[i, j]` gives the presence or absence of the influence of component `j` onto component `i`.

```julia
size(connec)        # Get the number of components.
```

# Remark
Since the connectiviy matrix is specified, the parameter `model.p` is not used.
"""
struct MarkovChainConnectivity
    model::MarkovChainModel
    θ::Matrix{Bool}
    function MarkovChainConnectivity(model::MarkovChainModel, θ::Matrix{Bool})
        if (size(θ)[1] != size(θ)[2])
            throw(DimensionMismatch("θ must be a square matrix"))
        end
        return new(model, θ)
    end
end

# constructors

# methods
Base.size(modelconnec::MarkovChainConnectivity) = size(modelconnec.θ)[1]

# print
function Base.show(io::IO, ::MIME"text/plain", modelconnec::MarkovChainConnectivity)
    return print(
        io,
        "μ = $(modelconnec.model.μ), λ = $(modelconnec.model.λ). The θ matrix is specified and has $(size(modelconnec)) nodes.",
    )
end

"""
    DiscreteTimeData(X)

A *binary discrete time data* of length `T` with `N` dimensions.

# Arguments
- `X::Matrix{Bool}`: `X[i, t]` gives the presence or absence of an event at time `t` for component `i`.

```julia
length(data)                # Get the time length of data
size(data)                  # Get the dimensions of data, i.e. `(N,T)`
data[:,range::UnitRange]    # Extract the data in the time interval `range`.
```
"""
struct DiscreteTimeData
    X::Matrix{Bool}
end

# constructors

# methods
Base.length(data::DiscreteTimeData) = size(data.X)[2]
Base.size(data::DiscreteTimeData) = size(data.X)
Base.getindex(data::DiscreteTimeData, range::UnitRange) = DiscreteTimeData(data.X[:, range])

# print
function Base.show(io::IO, ::MIME"text/plain", data::DiscreteTimeData)
    return print(
        io, "DiscreteTimeData with $(size(data)[1]) nodes and $(size(data)[2]) time steps."
    )
end

struct ErdosRenyiGraph
    N::Int64
    p::Float64
    function ErdosRenyiGraph(N::Int64, p::Float64)
        if (p < 0) || (p > 1)
            throw(ArgumentError("p is $p but must be in the range [0, 1]"))
        end
        return new(N, p)
    end
end

# constructors
ErdosRenyiGraph(N::Real, p::Real) = ErdosRenyiGraph(Int(N), Float64(p))

# methods

# print
function Base.show(io::IO, ::MIME"text/plain", model::ErdosRenyiGraph)
    return print(
        io,
        "Erdos Renyi graph with $(model.N) nodes and connexion parameter p = $(model.p).",
    )
end

## Auxiliary functions

"""
    mvw(μ, λ, p, r₊)

Compute the targets `m`, `v` and `w` corresponding to the parameters `μ`, `λ`, `p` and the ratio of excitatory components `r₊`. 

The mathematical expressions are:

```math
\\begin{aligned}
m &= \\frac{\\mu+(1-\\lambda)pr_-}{1-p(1-\\lambda)(r_+-r_-)},\\
v &= (1-\\lambda)^2 p(1-p)((m-r_-)^2+r_+r_-) , \\
w &= m(1-m)\\frac{1+4(1-\\lambda)^2p^2r_+r_-}{(1-p(1-\\lambda)(r_+-r_-))^2}.
\\end{aligned}
```

# Arguments
- `μ`: Spontaneous event probability 
- `λ`: Interaction parameter
- `p`: Edge probability
- `r₊`: ratio of excitatory components
"""
function mvw(μ, λ, p, r₊)
    if !(0 <= r₊ <= 1)
        throw(ArgumentError("r₊ is $r₊ but must be in the range [0, 1]"))
    end
    r₋ = 1 - r₊
    D = 1 - (1 - λ) * p * (r₊ - r₋)
    m = (μ + (1 - λ) * p * r₋) / D
    v = (1 - λ)^2 * p * (1 - p) * ((m - r₋)^2 + r₊ * r₋)
    w = m * (1 - m) * (1 + 4 * (1 - λ)^2 * p^2 * r₊ * r₋) / D^2
    return m, v, w
end

function mvw(model::MarkovChainModel, r₊)
    return mvw(model.μ, model.λ, model.p, r₊)
end

"""
    mvw_inf(modelconnec::MarkovChainConnectivity, excitatory::Vector{Bool})

Compute the limit of the three estimators ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` when ``T\\to \\infty``.

The mathematical expressions are:

```math
\\begin{aligned}
m_\\infty = \\overline{m^N}, \\
v_\\infty = \\sum_{i=1}^N (m^N_i - \\overline{m^N})^2 , \\
w_\\infty = \\frac{1}{N} \\sum_{i=1}^N (c^N_i)^2 m^N_i (1 - m^N_i) , 
\\end{aligned}
```
# Arguments
- `modelconnec::MarkovChainConnectivity`: a `MarkovChainModel` with a specified connectivity matrix `θ`.
- `excitatory::Vector{Bool}`: `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
"""
function mvw_inf(modelconnec::MarkovChainConnectivity, excitatory::Vector{Bool})
    N = length(excitatory)
    if N != size(modelconnec)
        throw(DimensionMismatch("size(modelconnec) must be equal to length(excitatory)"))
    end
    model = modelconnec.model
    θ = modelconnec.θ
    r₊ = sum(excitatory .== true) / N
    μ = model.μ
    λ = model.λ

    A = 1 / N * (θ .* transpose(-1 .+ 2 * excitatory))
    Q = inv(I - (1 - λ) * A)
    L⁻ = sum(A[:, excitatory .== false]; dims=2)
    mN = μ * Q * ones(N) - (1 - λ) * Q * L⁻
    c = sum(Q; dims=1)

    m_inf = mean(mN)
    v_inf = sum((mN .- m_inf) .^ 2)
    w_inf = mean(c .^ 2 .* (mN - mN .^ 2))
    return m_inf, v_inf, w_inf
end
function mvw_inf(modelconnec::MarkovChainConnectivity, N::Int, r₊::Float64)
    return mvw_inf(modelconnec, N2excitatory(N, r₊))
end

function N2excitatory(N::Int, r₊::Float64)
    if !(0 <= r₊ <= 1)
        throw(ArgumentError("r₊ is $r₊ but must be in the range [0, 1]"))
    end
    if (N < 1)
        throw(ArgumentError("N is $N but must be positive"))
    end

    N₊ = floor(Int, N * r₊)
    return [ones(Bool, N₊); zeros(Bool, N - N₊)]
end
