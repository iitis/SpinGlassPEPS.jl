using CUDA
using SpinGlassPEPS
using LinearAlgebra


using Test

my_tests = []
if CUDA.functional() && CUDA.has_cutensor()# && false
    push!(my_tests,
    "cuda/base.jl",
    "cuda/contractions.jl"
    )
end

push!(my_tests,
    "base.jl",
    "contractions.jl",
    "compressions.jl",
    "ising.jl"
)
for my_test in my_tests
    include(my_test)
end
