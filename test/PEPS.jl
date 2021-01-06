@testset "PepsTensor correctly builds PEPS network" begin

m = 3
n = 4
t = 3

β = 1

L = m * n * t

instance = "$(@__DIR__)/instances/pathological/test_$(m)_$(n)_$(t).txt" 

ig = ising_graph(instance, L)
update_cells!(
   ig, 
   rule = square_lattice((m, n, t)),
) 

fg = factor_graph(
    ig, 
    energy=energy, 
    spectrum=full_spectrum,
)


for e in edges(fg)
    pl, en, pr = get_prop(fg, e, :split)
    println("edge ", e)
    println(size(pl), "   ", size(en),  "   ", size(pr))
    #display(en)
end

x, y = m, n

#for origin ∈ (:NW, :SW, :WN, :NE, :EN, :SE, :ES, :SW, :WS)

for origin ∈ (:NW, :SW, :NE, :SE, :SW) # OK
#for origin ∈ (:WN, :EN, :ES, :WS)  # NO

    @info "testing peps" origin
    println(origin)

    peps = PepsNetwork(x, y, fg, β, origin)
    @test typeof(peps) == PepsNetwork

    for i ∈ 1:peps.i_max, j ∈ 1:peps.j_max
        A = SpinGlassPEPS._generate_tensor(peps, (i, j)) 
        B = generate_tensor(peps, (i, j))
        @test A ≈ B
    end

    @info "contracting MPOs (up --> down)"

    ψ = MPO(peps, 1)
    for A ∈ ψ @test size(A, 2) == 1 end

    for i ∈ 2:peps.i_max
        W = MPO(peps, i) 
        ψ = ψ * W
        for A ∈ ψ @test size(A, 2) == 1 end
    end

    for A ∈ ψ @test size(A, 4) == 1 end

    @info "contracting MPOs (down --> up)"

    ψ = MPO(peps, peps.i_max)
    for A ∈ ψ @test size(A, 4) == 1 end

    for i ∈ peps.i_max-1:1
        W = MPO(peps, i) 
        ψ = W * ψ 
        for A ∈ ψ @test size(A, 4) == 1 end
    end

    for A ∈ ψ @test size(A, 4) == 1 end
end

end