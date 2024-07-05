import Base.show, Base.length, Base.getindex, Base.size
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
end

# constructors

# methods

# print
show(io::IO, ::MIME"text/plain", model::MarkovChainModel) = print(io,"μ = $(model.μ), λ = $(model.λ). The θ matrix is not specified (but p = $(model.p)).")

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
end

# constructors

# methods
size(modelconnec::MarkovChainConnectivity) = size(modelconnec.θ)[1]

# print
show(io::IO, ::MIME"text/plain", modelconnec::MarkovChainConnectivity) = print(io,"μ = $(modelconnec.model.μ), λ = $(modelconnec.model.λ). The θ matrix is specified and has $(size(modelconnec)) nodes.")

"""
    DiscreteTimeData(X)

A *binary discrete time data* of length `T` with `N` dimensions.

# Arguments
- `X::Matrix{Bool}`: `X[i, t]` gives the presence or absence of an event at time `t` for component `i`.

```julia
legnth(data)                # Get the time length of data
size(data)                  # Get the dimensions of data, i.e. `(N,T)`
data[:,range::UnitRange]    # Extract the data in the time interval `range`.
```
"""
struct DiscreteTimeData
    X::Matrix{Bool}
end

# constructors

# methods
length(data::DiscreteTimeData) = size(data.X, 2)
size(data::DiscreteTimeData) = size(data.X)
getindex(data::DiscreteTimeData, range::UnitRange) = DiscreteTimeData(data.X[:,range])

# print
show(io::IO, ::MIME"text/plain", data::DiscreteTimeData) = print(io,"DiscreteTimeData with $(size(data)[1]) nodes and $(size(data)[2]) time steps.")


struct ErdosRenyiGraph
    N::Int64
    p::Float64
end

# constructors

# methods

# print
show(io::IO, ::MIME"text/plain", model::ErdosRenyiGraph) = print(io,"Erdos Renyi graph with $(model.N) nodes and connexion parameter p = $(model.p).")



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
    r₋ = 1 - r₊
    D = 1 - (1-λ)*p*(r₊ - r₋)
    m = (μ + (1-λ)*p*r₋) / D
    v = (1-λ)^2*p*(1-p)*( (m - r₋)^2 + r₊*r₋ )
    w = m*(1-m)*( 1 + 4*(1-λ)^2*p^2*r₊*r₋ ) / D^2
    return m,v,w
end

function mvw(model::MarkovChainModel, r₊)
    return mvw(model.μ, model.λ, model.p, r₊)
end

"""
    mvw_inf(modelconnec::MarkovChainConnectivity, Z::Vector{Bool})

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
- `Z::Vector{Bool}`: `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
"""
function mvw_inf(modelconnec::MarkovChainConnectivity, Z::Vector{Bool})
    model = modelconnec.model
    θ = modelconnec.θ
    N = length(Z)
    r₊ = sum(Z.==true)/N
    μ = model.μ
    λ = model.λ

    A = 1/N * (θ .* transpose(-1 .+ 2*Z))
    Q = inv(I - (1-λ)*A)
    L⁻ = sum(A[:,Z .== false], dims=2)
    mN = μ*Q*ones(N) - (1-λ)*Q*L⁻
    c = sum(Q, dims=1)

    m_inf = mean(mN)
    v_inf = sum((mN .- m_inf).^2)
    w_inf = mean(c.^2 .* (mN - mN.^2))
    return m_inf, v_inf, w_inf
end
mvw_inf(modelconnec::MarkovChainConnectivity, N::Int, r₊::Float64) = mvw_inf(modelconnec, N2Z(N,r₊))

function N2Z(N, r₊)
    N₊ = floor(Int,N*r₊)
    return [ones(Bool, N₊); zeros(Bool, N - N₊)]
end