mutable struct ConstraintGraph{N,T<:Integer}
    G::Graphs.SimpleDiGraph{T}
    constraints::Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}
end

function ConstraintGraph{N,T}() where {N,T<:Integer}
    G = Graphs.SimpleDiGraph{T}()
    constraints = Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}()
    return ConstraintGraph(G, constraints)
end

ConstraintGraph() = ConstraintGraph{3,Int}()

"""
    add_vertex(g::ConstraintGraph{N, T})

Add a new vertex to the constraint graph `g`.
"""
function add_vertex!(g::ConstraintGraph{N,T}) where {N,T<:Integer}
    Graphs.add_vertex!(g.G)
    return g
end

function add_constraint!(g::ConstraintGraph{N,T}, edge::Graphs.Edge{T}, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}
    Graphs.add_edge!(g.G, edge)
    (g.constraints[edge] = constraint)
    return g
end

"""
    add_constraint!(g::constraintGraph{N,T}, from::T, to::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}

Add a constraint between two vertices in the constraint graph `g`.
"""
add_constraint!(g::ConstraintGraph{N,T}, from::T, to::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer} = add_constraint!(g, Graphs.Edge(from, to), constraint)