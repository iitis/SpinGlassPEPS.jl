export Chimera, Lattice
export factor_graph, decompose_edges!
export Cluster, Spectrum
export rank_reveal

const SimpleEdge = LightGraphs.SimpleGraphs.SimpleEdge
const EdgeIter = Union{LightGraphs.SimpleGraphs.SimpleEdgeIter, Base.Iterators.Filter, Array}

mutable struct Cluster
    tag::Int
    vertices::Dict{Int,Int}
    edges::EdgeIter
    rank::Vector
    J::Matrix{<:Number}
    h::Vector{<:Number}

    function Cluster(ig::MetaDiGraph, v::Int, vertices::Dict, edges::EdgeIter)
        cl = new(v, vertices, edges)
        L = length(cl.vertices)

        cl.J = zeros(L, L)
        for e ∈ cl.edges
            i = cl.vertices[src(e)]
            j = cl.vertices[dst(e)] 
            cl.J[i, j] = get_prop(ig, e, :J)
        end

        rank = get_prop(ig, :rank)
        cl.rank = rank[1:L]

        cl.h = zeros(L)
        for (w, i) ∈ cl.vertices
            cl.h[i] = get_prop(ig, w, :h)
            cl.rank[i] = rank[w]
        end
        cl
    end
end

function MetaGraphs.filter_edges(ig::MetaGraph, v::Cluster, w::Cluster)
    edges = []
    for i ∈ keys(v.vertices)
        for j ∈ unique_neighbors(ig, i)
            if j ∈ keys(w.vertices)
                push!(edges, SimpleEdge(i, j))
            end
        end
    end
    edges
end

mutable struct Edge
    tag::NTuple
    edges::EdgeIter
    J::Matrix{<:Number}

    function Edge(ig::MetaDiGraph, v::Cluster, w::Cluster)
        ed = new((v.tag, w.tag))
        ed.edges = filter_edges(ig, v, w) 

        m = length(v.vertices)
        n = length(w.vertices)

        ed.J = zeros(m, n)
        for e ∈ ed.edges
            i = v.vertices[src(e)]
            j = w.vertices[dst(e)] 
            ed.J[i, j] = get_prop(ig, e, :J)
        end
        ed
    end
end

function factor_graph(m::Int, n::Int, hdir=:LR, vdir=:BT)
    dg = MetaDiGraph(SimpleDiGraph(m * n))
    set_prop!(dg, :order, (hdir, vdir))

    linear = LinearIndices((1:m, 1:n))
    for i ∈ 1:m
        for j ∈ 1:n-1
            v, w = linear[i, j], linear[i, j+1]
            hdir == :LR ? e = SimpleEdge(v, w) : e = SimpleEdge(w, v)
            add_edge!(dg, e)
            set_prop!(dg, e, :orientation, "horizontal")
        end
    end

    for i ∈ 1:n
        for j ∈ 1:m-1
            v, w = linear[j, i], linear[j+1, i]
            vdir == :BT ? e = SimpleEdge(v, w) : e = SimpleEdge(w, v)
            add_edge!(dg, e)
            set_prop!(dg, e, :orientation, "vertical")
        end
    end
    dg
end

function factor_graph(
    g::Model, 
    energy::Function=energy, 
    spectrum::Function=brute_force, 
    cluster::Function=unit_cell,
) 
    m, n, _ = g.size
    fg = factor_graph(m, n)

    for v ∈ vertices(fg)
        cl = cluster(g, v)
        set_prop!(fg, v, :cluster, cl)
        set_prop!(fg, v, :spectrum, spectrum(cl))
    end

    for e ∈ edges(fg)
        v = get_prop(fg, src(e), :cluster)
        w = get_prop(fg, dst(e), :cluster)

        edge = Edge(g.graph, v, w)
        set_prop!(fg, e, :edge, edge)
        set_prop!(fg, e, :energy, energy(fg, edge))
    end
    fg
end


function decompose_edges!(fg::MetaDiGraph, order=:PE, beta::Float64=1.0)
    set_prop!(fg, :tensorsOrder, order)

    for edge ∈ edges(fg)
        energies = get_prop(fg, edge, :energy)
        en, p = rank_reveal(energies)

        if order == :PE
            dec = (p, exp.(-beta .* en))
        else
            dec = (exp.(-beta .* en), p)
        end
        set_prop!(fg, edge, :decomposition, dec)
    end 
end
 

function rank_reveal(energy, order=:PE) # or E_then_P
    dim = order == :PE ? 1 : 2

    E = unique(energy, dims=dim)
    idx = indexin(eachslice(energy, dims=dim), collect(eachslice(E, dims=dim)))

    P = order == :PE ? zeros(size(energy, 1), size(E, 1)) : zeros(size(E, 2), size(energy, 2))
    
    for (i, elements) ∈ enumerate(eachslice(P, dims=dim))
        elements[idx[i]] = 1
    end

    order == :PE ? (P, E) : (E, P)
end 