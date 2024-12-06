## Main functions
"""
    estimators(data::DiscreteTimeData, Δ::Int=floor(Int,log(length(data))))

Compute the estimators ``\\hat{m}``, ``\\hat{v}`` and ``\\hat{w}`` on a data set `data` with tuning parameter `Δ`.

# Remark
If the value `Δ = 0` is chosen, then it is replaced by its default value `floor(Int,log(length(data)))`
"""
function estimators(
    data::DiscreteTimeData, Δ::Int=floor(Int, log(length(data)))
)::Tuple{Float64,Float64,Float64}
    if Δ == 0
        Δ = floor(Int, log(length(data)))
    end

    N, T = size(data)
    X = data.X

    Z_T = sum(X; dims=2)
    Z̄_T = mean(Z_T)
    m̂ = Z̄_T / T
    v̂ = (N) * (T + 1) * T^(-3) * (mean(Z_T .^ 2) - T / (T + 1) * (Z̄_T + Z̄_T^2))
    WΔ = 0
    for iter in 1:div(T, Δ)
        WΔ += (N) / T * (sum(X[:, (1 + (iter - 1) * Δ):(iter * Δ)]) / (N) - Δ * m̂)^2
    end
    W2Δ = 0
    for iter in 1:div(T, 2 * Δ)
        W2Δ +=
            (N) / T *
            (sum(X[:, (1 + (iter - 1) * 2 * Δ):(iter * 2 * Δ)]) / (N) - 2 * Δ * m̂)^2
    end
    ŵ = 2 * W2Δ - WΔ

    return (m̂, v̂, ŵ)
end

function estimators(
    data::DiscreteTimeData, Δvec::Vector{Int}
)::Tuple{Float64,Float64,Vector{Float64}}
    N, T = size(data)
    X = data.X

    Z_T = sum(X; dims=2)
    Z̄_T = mean(Z_T)
    m̂ = Z̄_T / T
    v̂ = (N) * (T + 1) * T^(-3) * (mean(Z_T .^ 2) - T / (T + 1) * (Z̄_T + Z̄_T^2))
    ŵ = Float64[]
    for Δ in Δvec
        if Δ == 0
            Δ = floor(Int, log(length(data)))
        end
        WΔ = 0
        for iter in 1:div(T, Δ)
            WΔ += (N) / T * (sum(X[:, (1 + (iter - 1) * Δ):(iter * Δ)]) / (N) - Δ * m̂)^2
        end
        W2Δ = 0
        for iter in 1:div(T, 2 * Δ)
            W2Δ +=
                (N) / T *
                (sum(X[:, (1 + (iter - 1) * 2 * Δ):(iter * 2 * Δ)]) / (N) - 2 * Δ * m̂)^2
        end
        push!(ŵ, 2 * W2Δ - WΔ)
    end

    return (m̂, v̂, ŵ)
end

"""
    fit(::Type{MarkovChainModel}, data::DiscreteTimeData, r₊; Δ::Int=floor(Int,log(length(data))))

Fit a `MarkovChainModel` to a data set `data` with asymptotic excitatory proportion `r₊`.

# Keyword argument
- `Δ`: tuning parameter of the method. Its default value is the floor of ``\\log(T)`` where ``T`` is the time length of `data`.
"""
function Distributions.fit(
    ::Type{MarkovChainModel},
    data::DiscreteTimeData,
    r₊::Float64;
    Δ::Int=floor(Int, log(length(data))),
)::MarkovChainModel
    m̂, v̂, ŵ = estimators(data, Δ)

    μ, λ, p = Φ(m̂, v̂, ŵ, r₊)
    return MarkovChainModel(μ, λ, p)
end

## Auxiliary functions
function ϕ(m::Float64, v::Float64, w_or_d::Float64, r₊::Float64)::Tuple{Float64,Float64}
    r₋ = 1 - r₊
    if abs(r₊ - r₋) < 1e-3
        ϕ₁ = w_or_d / (m * (1 - m)) - 1
    else
        ϕ₁ = (1 - w_or_d)^2 / (r₊ - r₋)^2
    end
    if ϕ₁ < 0
        @debug "Negative square ! Absolute value applied."
        ϕ₁ = abs(ϕ₁)
    end
    ϕ₂ = 1 + v / (((m - r₋)^2 + r₊ * r₋) * ϕ₁)
    return ϕ₁, ϕ₂
end

function Φ_aux(m::Float64, v::Float64, w_or_d::Float64, r₊::Float64)::Tuple{Float64,Float64,Float64}
    r₋ = 1 - r₊
    ϕ₁, ϕ₂ = ϕ(m, v, w_or_d, r₊)
    Φ₁ = m * (1 - (r₊ - r₋) * sqrt(ϕ₁)) - r₋ * sqrt(ϕ₁)
    if ϕ₂ == Inf
        Φ₂ = 1.0
    else
        Φ₂ = 1 - ϕ₂ * sqrt(ϕ₁)
    end
    Φ₃ = 1 / ϕ₂
    distance2admissibleset(Φ₁, Φ₂, Φ₃) == 0 ||
        @debug "Estimators were clipped to admissible values."
    return projection2admissibleset(Φ₁, Φ₂, Φ₃)
end

function Φ(m::Float64, v::Float64, w::Float64, r₊::Float64)::Tuple{Float64,Float64,Float64}
    r₋ = 1 - r₊
    if abs(r₊ - r₋) < 1e-3
        return Φ_aux(m, v, w, r₊)
    end

    κ = (r₊ - r₋)^2 * w / (m * (1 - m))

    if abs(κ - 4 * r₊ * r₋) < 1e-3
        d = (8 * r₊ * r₋)^(-1)
        return Φ_aux(m, v, d, r₊)
    end

    discriminant = κ - 4 * r₊ * r₋ + (4 * r₊ * r₋)^2
    if discriminant <= 0
        @debug "Negative discriminant ! Clipped to 0."
        d = 4 * r₊ * r₋ / (4 * r₊ * r₋ - κ)
        return Φ_aux(m, v, d, r₊)
    end

    d₊ = (4 * r₊ * r₋ + sqrt(discriminant)) / (4 * r₊ * r₋ - κ)
    d₋ = (4 * r₊ * r₋ - sqrt(discriminant)) / (4 * r₊ * r₋ - κ)
    if r₊ > 0.5
        d = d₋
        0 < d < 1 || @debug "Careful. r₊=" *
            string(r₊) *
            ", d₋=" *
            string(d₋) *
            " and d₊=" *
            string(d₊)
        return Φ_aux(m, v, d, r₊)
    end
    if 4 * r₊ * r₋ - κ < 0
        d = d₋
        d > 1 || @debug "Careful. 4*r₊*r₋ - κ=" *
            string(4 * r₊ * r₋ - κ) *
            ", d₋=" *
            string(d₋) *
            " and d₊=" *
            string(d₊)
        return Φ_aux(m, v, d, r₊)
    end
    if d₊ > 2 * r₋
        d = d₋
        d > 1 || @debug "Careful. d₊ > 2*r₋, d₋=" * string(d₋) * " and d₊=" * string(d₊)
        return Φ_aux(m, v, d, r₊)
    end

    μ₊, λ₊, p₊ = Φ_aux(m, v, d₊, r₊)
    distfix₊ = sum((mvw(μ₊, λ₊, p₊, r₊) .- (m, v, w)) .^ 2)
    μ₋, λ₋, p₋ = Φ_aux(m, v, d₋, r₊)
    distfix₋ = sum((mvw(μ₋, λ₋, p₋, r₊) .- (m, v, w)) .^ 2)
    whichestim = argmin([distfix₊, distfix₋])
    @debug "Identifiability problem. The triplet with the lowest fixed point distance is returned."
    return (whichestim == 1) .* (μ₊, λ₊, p₊) .+ (whichestim == 2) .* (μ₋, λ₋, p₋)
end

function distance2admissibleset(μ::Float64, λ::Float64, p::Float64)::Float64
    d1 = -λ * (λ < 0) + (λ - 1) * (λ > 1)
    d2 = -p * (p < 0) + (p - 1) * (p > 1)
    d3 = -μ * (μ < 0) + (μ - λ) * (μ > λ)
    return d1 + d2 + d3
end

function projection2admissibleset(μ::Float64, λ::Float64, p::Float64)::Tuple{Float64,Float64,Float64}
    λ = min(1, max(0, λ))
    p = min(1, max(0, p))
    μ = min(λ, max(0, μ))
    return (μ, λ, p)
end
