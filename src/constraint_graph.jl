mutable struct ConstraintGraph{N,T<:Integer}
    G::Graphs.SimpleDiGraph{T}
    constraints::Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}
    filtered::Bool # flag that indicates if the constraints have been filtered

    function ConstraintGraph(g::Graphs.SimpleDiGraph{T}, constraints::Dict{Graphs.Edge{T},DisplacementGroupComposition{N}}) where {N,T<:Integer}
        new{N,T}(g, constraints, false)
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
    g.filtered = false
    return g
end

"""
    add_constraint!(g::ConstraintGraph{N,T}, src::T, dst::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}

Add a constraint between two vertices in the constraint graph `g`. Return true if the constraint was added successfully, otherwise return false.
"""
function add_constraint!(g::ConstraintGraph{N,T}, src::T, dst::T, constraint::DisplacementGroupComposition{N}) where {N,T<:Integer}
    if Graphs.add_edge!(g.G, src, dst)
        g.constraints[Graphs.Edge(src, dst)] = constraint
        g.filtered = false
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
"""
    filter_constraints!(g::ConstraintGraph{N,T}) where {N,T<:Integer}

Filter the constraints of the constraint graph g along the basis cycles of g. See Thomas 1991, section 4.4.3
"""
function filter_constraints!(g::ConstraintGraph{N,T}) where {N,T<:Integer}
    if g.filtered
        return g
    end

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

            if length(updated_constraints) == 0
                repeat = false
            end

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

    g.filtered = true

    return g
end

filter_constrints(g::ConstraintGraph{N,T}) where {N,T<:Integer} = filter_constraints!(deepcopy(g))


"""
    all_equiv_constraints(g::ConstraintGraph, src::Integer)

Compute all equivalent constraints of the constraint graph g between the vertex src and all other vertices in the graph. 
This is achieved by filtering the constraints in the constraint graph and then composing the constraints along shortest paths from src to the other vertices.
The shortest paths are calculated using Dijkstra's algorithm.
"""
function all_equiv_constraints!(g::ConstraintGraph{N,T}, src::Integer) where {N,T<:Integer}
    # filter constraints
    filter_constraints!(g)

    f = deepcopy(g)

    # we are interested in paths that don't care about the direction of the edges -> add reversed edges with inverse constraints
    for e in Graphs.edges(g.G)
        add_constraint!(f, e.dst, e.src, inv(g.constraints[e]))
    end

    # compute the shortest paths
    dijk = Graphs.dijkstra_shortest_paths(f.G, src)

    # compute the equivalent constraints by following the paths computed by Dijkstra's algorithm
    # this can be done in O(n) time 
    equiv_constraints = [DisplacementGroupComposition([IdentityGroup(N)]) for _ in 1:Graphs.nv(f.G)]
    for v in sort(collect(1:Graphs.nv(f.G)), by=v -> dijk.dists[v])
        if dijk.dists[v] == 0
            continue
        end

        # equivalent constraint is the composition of the equivalent constraint of the parent vertex of v along a shortest path with the filtered constraint of the edge connecting the parent vertex with v
        equiv_constraints[v] = equiv_constraints[dijk.parents[v]] * f.constraints[Graphs.Edge(dijk.parents[v], v)]
    end

    return equiv_constraints
end

# function that checks if a kinematic constraint graph is interlocking when a given frame is fixed.
# This is done by adding a global reference frame (a new vertex) and fixing the frame vertices relative to the global reference frame (by imposing constraints = I).
# Then the equivalent constraints between the global reference frame and all other vertices are computed.
# If all equivalent constraints are trivial, the constraint graph is interlocking. (The allowed motions of the blocks are a subset of the constraints).
function is_interlocking(g::ConstraintGraph{N,T}, frame::AbstractVector{<:Integer}) where {N,T<:Integer}
    if !issubset(frame, collect(1:Graphs.nv(g.G)))
        throw(ArgumentError("Frame vertices are not a subset of the vertices of the constraint graph"))
    end

    f = deepcopy(g)

    add_vertex!(f) # add vertex representing global reference frame
    for v in frame
        add_constraint!(f, Graphs.nv(f.G), v, DisplacementGroupComposition([IdentityGroup(N)])) # fix frame vertices relative to global reference frame
    end

    # check if the constraint graph is interlocking
    equiv_constraints = all_equiv_constraints!(f, Graphs.nv(f.G))
    return all(constraint -> is_trivial!(constraint) && constraint.factors[1] == IdentityGroup(N), equiv_constraints)
end

all_equiv_constraints(g::ConstraintGraph{N,T}, src::Integer) where {N,T<:Integer} = all_equiv_constraints!(deepcopy(g), src)

#################### plotting of constraint graphs using GraphMakie
function GraphMakie.graphplot(g::ConstraintGraph, kwargs...)
    G = g.G
    labels = repr.(g.constraints[e] for e in Graphs.edges(G))
    f, ax, p = graphplot(G, elabels=labels, kwargs...)

    return f, ax, p
end