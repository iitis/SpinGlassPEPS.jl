@testset "cuda" begin
D = 10
Dcut = 5
d = 2
sites = 5

T = Float16
@testset "creation of random MPS on cuda" begin
    ψ = randn(MPS{T}, sites, D, d)

    @test length(ψ) == sites
    @test ψ == ψ
    @test ψ ≈ ψ
    @test eltype(ψ) == T
    @test typeof(ψ.tensors[1]) <: CuArray
    @test copy(ψ) == ψ
end


@testset "random MPO on cuda" begin
    H = randn(MPO{T}, sites, D, d)
    @test typeof(H.tensors[1]) <: CuArray
end

@testset "contraction on cuda" begin
    ψ = randn(MPS{T}, sites, D, d)
    H = randn(MPO{T}, sites, D, d)

    @test typeof((H*ψ).tensors[1]) <: CuArray
end
end