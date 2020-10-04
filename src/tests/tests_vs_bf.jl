using NPZ
using Plots
using Test

include("../notation.jl")
include("../peps.jl")

file = "examples.npz"
β = 2.
j = 10
examples = 100
r=2

if false
    file = "examples2.npz"
    β = 2.
    j = 10
    examples = 100
    r = 1
end

if false
file = "energies_and_matrix_only.npz"
    β = 2.
    j = 25
    examples = 1
    r = 2
end

T = Float64

data = npzread("./data/"*file)
println(file)

function v2energy(M::Matrix{T}, v::Vector{Int}) where T <: AbstractFloat
    d =  diag(M)
    M = M .- diagm(d)

    transpose(v)*M*v + transpose(v)*d
end


for k in 1:examples
    println("SAMPLE = ", k)
    QM = data["Js"][k,:,:]

    states = 0
    try
        states = data["states"][k,:,:]
    catch
        0
    end
    ens = data["energies"][k,:,:]

    function M2Qubbo_els(M::Matrix{Float64}, T::Type = Float64)
        qubo = Qubo_el{T}[]
        s = size(M)
        for i in 1:s[1]
            for j in i:s[2]
                if (M[i,j] != 0.) | (i == j)
                    x = T(M[i,j])
                    q = Qubo_el{T}((i,j), x)
                    push!(qubo, q)
                end
            end
        end
        qubo
    end

    qubo = M2Qubbo_els(QM, T)

    grid = [1 2 3 4 5; 6 7 8 9 10; 11 12 13 14 15; 16 17 18 19 20; 21 22 23 24 25]


    ses = solve(qubo, grid, r*j; β = T(β), χ = 0, threshold = T(0.))

    count = copy(j)
    for i in 1:j

        v = ses[i].spins

        if !(v2energy(QM, v) ≈ -ens[i])
            println("exact")
            println("n.o. state = ", i)
            println("energies (peps, bf)", (v2energy(QM, v), -ens[i]))
            count = count - 1

            try
                v1 = Int.(states[i,:])
                println(v1)
                println(v)
            catch
                0
            end
        end
    end

    if count != j
        println("exact n.o. corresponding energies $(count), should be $j")
    end

    probs = [ses[i].objective for i in 1:j]
    if k == 1
        energies_exact = [ens[i] for i in 1:j]

        ps = exp.(-energies_exact*β)

        y = [v2energy(QM, ses[i].spins) for i in 1:j]
        y = exp.(y*β)
        p_theoretical = ps./sum(ps)
        A = [p_theoretical, y./sum(y), probs]

        plot(A, label = ["bf" "peps M" "peps"], yaxis = :log)
        savefig("./pics/$(file)_$(k)_myplot.pdf")
    end


    χ = 2
    ses_a = solve(qubo, grid, r*j; β = T(β), χ = χ, threshold = T(1e-10))


    count_a = copy(j)
    for i in 1:j
        v = ses_a[i].spins

        if !(v2energy(QM, v) ≈ -ens[i])
            println("approx")
            println("n.o. state = ", i)
            println("energies (peps,bf)", (v2energy(QM, v), -ens[i]))
            count_a = count_a - 1
            try
                v1 = Int.(states[i,:])
                println(v1)
                println(v)
            catch
                0
            end
        end
    end

    if count_a != j
        println("approx n.o. corresponding energies = $(count_a), should be $j")
    end

    probs_a = [ses_a[i].objective for i in 1:j]

    println("max diff exact and approx = ", maximum(abs.(probs .- probs_a)))
end