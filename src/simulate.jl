import Base.rand

## Main functions
function rand(graph::ErdosRenyiGraph)::Matrix{Bool}
    N = graph.N
    p = graph.p
    Adj = rand(Bernoulli(p), (N, N))
    return(Adj)
end

"""
    rand(modelconnec::MarkovChainConnectivity, Z::Vector{Bool}, T::Int)

Simulate a realization of a *Markov Chain model* with a specified connectivity matrix.

# Arguments
- `modelconnec::MarkovChainConnectivity`: a `MarkovChainModel` with a specified connectivity matrix `θ`.
- `Z::Vector{Bool}`: `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
- `T::Int`: Time length of the simulation.

```jldoctest
using MeanFieldGraph
model = MarkovChainModel(.5, .5, .5)
θ = [[1 1];[0 1]]
modelconnec = MarkovChainConnectivity(model,θ)
Z = [true, false]

using Random
Random.seed!(1)
data = rand(modelconnec, Z, 10)
data.X

# output

2×10 Matrix{Bool}:
 1  0  0  0  1  1  1  1  1  1
 1  0  1  1  0  1  0  1  1  0
```
"""
function rand(modelconnec::MarkovChainConnectivity, Z::Vector{Bool}, T::Int)::DiscreteTimeData
    N = length(Z)
    output = Matrix{Bool}(undef, N, T)
    current_value = stationary_initial_condition(modelconnec, Z)
    output[:,1] = current_value
    for t in 1:(T-1)
        forward_simulation!(current_value, modelconnec, Z)
        output[:,t+1] = current_value
    end
    return DiscreteTimeData(output)
end
rand(modelconnec::MarkovChainConnectivity, N::Int, r₊::Float64, T::Int)::DiscreteTimeData = rand(modelconnec, N2Z(N, r₊), T)


"""
    rand(model::MarkovChainModel, Z::Vector{Bool}, T::Int)

Simulate a realization of a *Markov Chain model* without a specified connectivity matrix.

# Arguments
- `model::MarkovChainModel`: contains the parameters `μ`, `λ` and `p`.
- `Z::Vector{Bool}`: `true` coordinates correspond to excitatory components and `false` coordinates correspond to inhibitory components.
- `T::Int`: Time length of the simulation.

# Remark
It generates a connectivity matrix according to an Erdos-Rényi graph with parameter `p` and then apply the method `rand(modelconnec::MarkovChainConnectivity, Z::Vector{Bool}, T::Int)`.
```
"""
function rand(model::MarkovChainModel, Z::Vector{Bool}, T::Int)::DiscreteTimeData
    graph = ErdosRenyiGraph(length(Z), model.p)
    θ = rand(graph)
    rand(MarkovChainConnectivity(model,θ), Z, T)
end
rand(model::MarkovChainModel, N::Int, r₊::Float64, T::Int)::DiscreteTimeData = rand(model, N2Z(N, r₊), T)


## Auxiliary functions
function stationary_initial_condition(modelconnec::MarkovChainConnectivity, Z::Vector{Bool})::Vector{Bool}
    history_values = Dict{Int, Bool}[]
    history_child = Vector{Int}[]
    history_parent = Vector{Int}[]
    unknown_nodes = collect(1:size(modelconnec))
    while !isempty(unknown_nodes)
        values, child, parent = backward_step(unknown_nodes, modelconnec)
        push!(history_values, values)
        push!(history_child, child)
        push!(history_parent, parent)
        unknown_nodes = unique(parent)
    end
    pop!(history_child); pop!(history_parent);
    parent_values = pop!(history_values)

    if isempty(history_child)
        return dicttovec(parent_values)
    end

    for _ in 1:length(history_child)
        child = pop!(history_child)
        parent = pop!(history_parent)
        child_values = pop!(history_values)
        for (index, i) in enumerate(child)
            j = parent[index]
            if Z[j] == true
                child_values[i] = parent_values[j]
            else
                child_values[i] = 1-parent_values[j]
            end
        end
        parent_values = child_values
    end
    return dicttovec(parent_values)
end

function backward_step(unknown_nodes::Vector{Int}, modelconnec::MarkovChainConnectivity)::Tuple{Dict{Int, Bool}, Vector{Int}, Vector{Int}}
    model = modelconnec.model
    N = size(modelconnec)
    λ = model.λ
    β = model.μ / λ
    θ = modelconnec.θ
    values = Dict{Int, Bool}()
    remaining_nodes = Int[]
    parent_nodes = Int[]
    for i in unknown_nodes
        if rand(Bernoulli(λ))
            values[i] = rand(Bernoulli(β))
        else
            j = rand(DiscreteUniform(1,N))
            if θ[i, j] == 0
                values[i] = 0
            else
                push!(parent_nodes, j)
                push!(remaining_nodes, i)
            end
        end
    end
    return values, remaining_nodes, parent_nodes
end

function forward_simulation!(current_value::Vector{Bool}, modelconnec::MarkovChainConnectivity, Z::Vector{Bool})
    past_value = copy(current_value)

    model = modelconnec.model
    N = size(modelconnec)
    λ = model.λ
    β = model.μ / λ
    θ = modelconnec.θ
    for i in 1:N
        if rand(Bernoulli(λ))
            current_value[i] = rand(Bernoulli(β))
        else
            j = rand(DiscreteUniform(1,N))
            if θ[i, j] == 0
                current_value[i] = 0
            else
                if Z[j] == true
                    current_value[i] = past_value[j]
                else
                    current_value[i] = 1-past_value[j]
                end
            end
        end
    end
end


function dicttovec(dict::Dict{Int, Bool})::Vector{Bool}
    n = length(dict)
    output = Vector{Bool}(undef, n)
    for index in 1:n
        output[index] = dict[index]
    end
    return output
end