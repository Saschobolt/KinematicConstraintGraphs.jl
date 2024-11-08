mutable struct ConstraintGraph{N,T<:Integer}
    G::Graphs.SimpleDiGraph{T}
    constraints::Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}

    function ConstraintGraph(g::Graphs.SimpleDiGraph{T}, constraints::Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}) where {N,T<:Integer}
        new{N,T}(g, constraints)
    end
end

function ConstraintGraph{N,T}() where {N,T<:Integer}
    G = Graphs.SimpleDiGraph{T}()
    constraints = Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}()
    return ConstraintGraph(G, constraints)
end

ConstraintGraph() = ConstraintGraph{3,Int}()

Base.show(io::IO, g::ConstraintGraph{N,T}) where {N,T} = print(io, "{$(Graphs.nv(g.G)), $(Graphs.ne(g.G))} $T kinematic constraint graph in dimension $N")

function Base.getindex(g::ConstraintGraph{N,T}, inds...) where {N,T<:Integer}
    F = getindex(g.G, inds...)
    constraints = Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}()
    for e in Graphs.edges(F)
        constraints[e] = g.constraints[e]
    end
    return ConstraintGraph(F, constraints)
end

"""
    add_vertex(g::ConstraintGraph{N, T})

Add a new vertex to the constraint graph `g`.
"""
function add_vertex!(g::ConstraintGraph{N,T}) where {N,T<:Integer}
    Graphs.add_vertex!(g.G)
    return g
end

"""
    add_constraint!(g::ConstraintGraph{N,T}, src::T, dst::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}

Add a constraint between two vertices in the constraint graph `g`. Return true if the constraint was added successfully, otherwise return false.
"""
function add_constraint!(g::ConstraintGraph{N,T}, src::T, dst::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}
    if Graphs.add_edge!(g.G, src, dst)
        g.constraints[Graphs.Edge(src, dst)] = constraint
        return true
    end

    return false
end

add_constraint!(g::ConstraintGraph{N,T}, src::T, dst::T, constraint::AbstractDisplacementGroup{N}) where {N,T<:Integer} = add_constraint!(g, src, dst, DisplacementGroupComposition([constraint]))


######### blocks of a graph
"""
    blocks(g::Union{Graphs.SimpleGraph,Graphs.SimpleDiGraph})

Compute the blocks (2-connected components) of a graph g and return them as a vector of vectors of vertices.
"""
function blocks(g::Union{Graphs.SimpleGraph,Graphs.SimpleDiGraph}) # TODO: improve efficiency of implementation. See maybe https://en.wikipedia.org/wiki/Biconnected_component
    cutpoints = Graphs.articulation(Graphs.SimpleGraph(g))
    # connected components of graph with edges of cutpoints removed
    f = copy(g)
    rem_edges = filter(e -> any(x -> x in [e.src, e.dst], cutpoints), collect(Graphs.edges(f)))
    for e in rem_edges
        Graphs.rem_edge!(f, e)
    end
    comps = setdiff(Graphs.connected_components(f), [[v] for v in cutpoints])

    # find blocks by checking which cutpoints have connection to which connected components
    for component in comps
        cutpoints_in_component = filter(x -> length(intersect(Graphs.all_neighbors(g, x), component)) > 0, cutpoints)
        append!(component, cutpoints_in_component)
        sort!(component)
    end

    return sort(comps)
end

####### filter constraints in kinematic constraint graph as in Thomas 1991 4.4.3
function filter_constraints!(g::ConstraintGraph{N,T}) where {N,T<:Integer}
    B = blocks(g.G)
    filter!(b -> length(b) > 1, B) # remove bridges from the constraint graph as those constraints do not propagate through the graph

    for block in B
        repeat = true
        while repeat
            f = g[block]
            cycles = Graphs.cycle_basis(Graphs.Graph(f.G))

            # dict that contains the propagation of the constraints through the graph. propagation_dict[(i, e)] is the result of filtering the constraint of the edge e through the cycle i
            propagation_dict = Dict{Tuple{Int,Graphs.Edge{T}},DisplacementGroupComposition{N}}()

            # check node consistency
            for (i, C) in enumerate(cycles)
                for n in eachindex(C)
                    # nth constraint in the cycle C
                    src = block[C[n]] # vertex of g that corresponds to C[n]
                    dst = block[C[mod1(n + 1, length(C))]] # next vertex in the cycle
                    if Graphs.Edge(src, dst) in Graphs.edges(g.G)
                        e = Graphs.Edge(src, dst)
                        direction = 1 # direction of the edge e in the cycle is forwards
                    else
                        e = Graphs.Edge(dst, src)
                        direction = -1 # direction of the edge e in the cycle is backwards
                    end

                    # filter the constraint of the edge e through the cycle C
                    comp = DisplacementGroupComposition([IdentityGroup(3)]) # composition that is intersected with the constraint of e
                    for m in vcat(n+1:length(C)..., 1:n-1...)
                        src2 = block[C[m]] # vertex of g that corresponds to C[m]
                        dst2 = block[C[mod1(m + 1, length(C))]] # next vertex in the cycle
                        if Graphs.Edge(src2, dst2) in Graphs.edges(g.G)
                            # direction of the edge e2 in the cycle is forwards
                            e2 = Graphs.Edge(src2, dst2)
                            comp = inv(g.constraints[e2]) * comp # constraint is multiplied inversely
                        else
                            # direction of the edge e2 in the cycle is backwards
                            e2 = Graphs.Edge(dst2, src2)
                            comp = g.constraints[e2] * comp # constraint is multiplied
                        end
                    end

                    if is_trivial!(comp)
                        if direction == 1
                            propagation_dict[(i, e)] = intersection(g.constraints[e], comp)
                        elseif direction == -1
                            propagation_dict[(i, e)] = intersection(g.constraints[e], inv(comp))
                        end
                    end
                end
            end

            # check arc consistency
            key_vec = collect(keys(propagation_dict))
            updated_constraints = unique([x[2] for x in key_vec])
            for e in updated_constraints
                # intersect constraint over all cycles it is a part of
                sol = DisplacementGroupComposition([SpecialEuclideanGroup(N)])
                for i in map(x -> x[1], filter(key -> key[2] == e, key_vec))
                    if !is_trivial!(propagation_dict[(i, e)])
                        continue
                    end
                    sol = intersection(sol, propagation_dict[(i, e)])
                end

                if sol != g.constraints[e]
                    g.constraints[e] = sol
                    repeat = true
                else
                    repeat = false
                end
            end
        end
    end

    return g
end