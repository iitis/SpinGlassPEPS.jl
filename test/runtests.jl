using CUDA
using SpinGlassPEPS
using LinearAlgebra
using TensorOperations
using TensorCast
using LightGraphs
using MetaGraphs
using Random
using Logging
using Statistics
using NPZ

disable_logging(LogLevel(1))

using Test

include("test_helpers.jl")

my_tests = []
if CUDA.functional() && CUDA.has_cutensor() && false
    CUDA.allowscalar(false)
    include("cuda/cu_test_helpers.jl")
    push!(my_tests,
    "cuda/base.jl",
    "cuda/contractions.jl",
    "cuda/compressions.jl",
    "cuda/spectrum.jl"
    )
end

include("test_helpers.jl")
push!(my_tests,
"MPS_search.jl",
#=
    "base.jl",
    "utils.jl",
    "contractions.jl",
    "compressions.jl",
    "identities.jl",
    "ising.jl",
    "indexing.jl",
    "searchMPS.jl",
    "MPS_search.jl",
    "factor.jl",
    "PEPS.jl",
    "contract.jl",
    "indexing.jl",
    =#
)

for my_test in my_tests
    include(my_test)
end
