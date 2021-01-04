export NetworkGraph, PepsNetwork
export generate_tensor, MPO

mutable struct NetworkGraph
    factor_graph::MetaDiGraph
    nbrs::Dict
    β::Number

    function NetworkGraph(factor_graph::MetaDiGraph, nbrs::Dict, β::Number)
        ng = new(factor_graph, nbrs, β) 

        count = 0
        for v ∈ vertices(ng.factor_graph), w ∈ ng.nbrs[v]
            if has_edge(ng.factor_graph, v, w) count += 1 end
        end

        mc = ne(ng.factor_graph)
        if count < mc
            error("Error: $(count) < $(mc)") 
        end
        ng
    end
end

# This still does not work (mixing @eval and @cast gives unexpected results)
function generate_tensor_general(ng::NetworkGraph, v::Int)
    loc_en = get_prop(ng.factor_graph, v, :loc_en)
    tensor = exp.(-ng.β .* loc_en)

    for w ∈ ng.nbrs[v]
        n = max(1, ndims(tensor)-1)
        s = :(@ntuple $n i)

        if has_edge(ng.factor_graph, w, v)
            pw, e, pv = get_prop(ng.factor_graph, w, v, :decomposition)
            @eval @cast tensor[σ, s, γ] |= tensor[σ, s] * pv[γ, σ]

        elseif has_edge(ng.factor_graph, v, w)
            pv, e, pw = get_prop(ng.factor_graph, v, w, :decomposition)
            @eval @cast tensor[σ, s, γ] |= tensor[σ, s] * pv[σ, γ]
        else 
            pv = ones(size(loc_en)...)
            @eval @cast tensor[σ, s, γ] |= tensor[σ, s] * pv[σ, γ]
        end
    end
    tensor
end

function generate_tensor(ng::NetworkGraph, v::Int)
    loc_en = get_prop(ng.factor_graph, v, :loc_en)
    ten_loc = exp.(-ng.β .* loc_en)

    p_list = Dict()
    for (i,w) ∈ enumerate(ng.nbrs[v])    
        if has_edge(ng.factor_graph, w, v)
            pw, e, pv = get_prop(ng.factor_graph, w, v, :decomposition)
            pv = pv'
        elseif has_edge(ng.factor_graph, v, w)
            pv, e, pw = get_prop(ng.factor_graph, v, w, :decomposition)
        else 
            pv = ones(length(loc_en), 1)
        end
        push!(p_list, i => pv)
    end
    @cast tensor[l, u, r, d, σ] |= ten_loc[σ] * p_list[1][σ, l] * p_list[2][σ, u] * p_list[3][σ, r] * p_list[4][σ, d] 
    tensor
end

function generate_tensor(ng::NetworkGraph, v::Int, w::Int)
    if has_edge(ng.factor_graph, w, v)
        _, e, _ = get_prop(ng.factor_graph, w, v, :decomposition)
        tensor = exp.(-ng.β .* e') 
    elseif has_edge(ng.factor_graph, v, w)
        _, e, _ = get_prop(ng.factor_graph, v, w, :decomposition)
        tensor = exp.(-ng.β .* e) 
    else 
        tensor = ones(1, 1)
    end
    tensor
end

mutable struct PepsNetwork
    size::NTuple{2, Int}
    map::Dict
    network_graph::NetworkGraph
    orientation::Symbol
    i_max::Int
    j_max::Int

    function PepsNetwork(m::Int, n::Int, fg::MetaDiGraph, β::Number, orientation::Symbol)
        pn = new((m, n)) 
        pn.map, pn.i_max, pn.j_max = LinearIndices(m, n, orientation)

        nbrs = Dict()
        for i ∈ 1:m, j ∈ 1:n
            push!(nbrs, 
            pn.map[i, j] => (pn.map[i, j-1], pn.map[i-1, j], 
                             pn.map[i, j+1], pn.map[i+1, j]))
        end
        pn.network_graph = NetworkGraph(fg, nbrs, β)
        pn
    end
end

generate_tensor(pn::PepsNetwork, m::NTuple{2,Int}) = generate_tensor(pn.network_graph, pn.map[m])
generate_tensor(pn::PepsNetwork, m::NTuple{2,Int}, n::NTuple{2,Int}) = generate_tensor(pn.network_graph, pn.map[m], pn.map[n])

function MPO(peps::PepsNetwork, i::Int; type::DataType=Float64)
    n = peps.j_max
    ψ = MPO(type, n)
    
    for j ∈ 1:n
        A = generate_tensor(peps, (i, j))
        @reduce B[l, u, r ,d] |= sum(σ) A[l, u, r, d, σ]
        ψ[j] = B
    end

    for j ∈ 1:n-1
        ten = generate_tensor(peps, (i, j), (i, j+1))
        A = ψ[j]
        @tensor B[l, u, r, d] := A[l, u, r̃, d] * ten[r̃, r]
        ψ[j] = B
    end
    ψ
end