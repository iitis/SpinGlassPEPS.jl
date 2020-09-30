include("peps.jl")
include("notation.jl")

function initialize_mps(l::Int, physical_dims::Int =  2, T::Type = Float64)
    [ones(T, 1,1,physical_dims) for _ in 1:l]
end


function make_ones(T::Type = Float64)
    d = 2
    ret = zeros(T, 1,1,d,d)
    for j in 1:d
        ret[1,1,j,j] = T(1.)
    end
    ret
end

function make_ones_inside(T::Type = Float64)
    d = 2
    ret = zeros(T, d,d,d,d)
    for i in 1:d
        for j in 1:d
            ret[i,i,j,j] = T(1.)
        end
    end
    ret
end

function T_with_B(l::Bool = false, r::Bool = false, T::Type = Float64)
    d = 2
    ret = zeros(T, d,d,d,d)
    for i in 1:d
        for j in 1:d
            for k in 1:d
                for l in 1:d
                    ret[i,j,k,l] = T(i==j)*T(k==l)*T(j==l)
                end
            end
        end
    end
    if l
        return sum(ret, dims = 1)
    elseif r
        return sum(ret, dims = 2)
    end
    return ret
end

function T_with_C(Jb::T, l::Bool = false, r::Bool = false) where T <: AbstractFloat
    d = 2
    ret = zeros(T, d,d,d,d)
    for i in 1:d
        for j in 1:d
            for k in 1:d
                for l in 1:d
                    ret[i,j,k,l] = T(i==j)*T(k==l)*exp(Jb*ind2spin(i)*ind2spin(k))
                end
            end
        end
    end
    if l
        return sum(ret, dims = 1)
    elseif r
        return sum(ret, dims = 2)
    end
    return ret
end

function add_MPO!(mpo::Vector{Array{T, 4}}, i::Int, i_n::Vector{Int}, qubo::Vector{Qubo_el{T}}, β::T) where T<: AbstractFloat
    k = minimum([i, i_n...])
    l = maximum([i, i_n...])
    for j in k:l
        mpo[j] = make_ones_inside(T)
    end
    mpo[i] = T_with_B(i==k, i==l, T)
    for j in i_n
        J = JfromQubo_el(qubo, i,j)
        mpo[j] = T_with_C(J*β, j==k, j==l)
    end
    mpo
end

function add_phase!(mps::Vector{Array{T, 3}}, qubo::Vector{Qubo_el{T}}, β::T) where T<: AbstractFloat
    d = size(mps[1], 3)
    for i in 1:length(mps)
        # we have a twice from the scalar product
        h = JfromQubo_el(qubo, i,i)/2
        for j in 1:d
            mps[i][:,:,j] = mps[i][:,:,j]*exp(ind2spin(j)*β*h)
        end
    end
end
