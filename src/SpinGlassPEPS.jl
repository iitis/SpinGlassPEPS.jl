module SpinGlassPEPS
    using LinearAlgebra
    using Requires
    using TensorOperations, TensorCast
    using LowRankApprox
    using LightGraphs
    using MetaGraphs
    using CSV

    include("base.jl")
    include("compressions.jl")
    include("contractions.jl")   
    include("ising.jl")
    include("graph.jl")
    include("PEPS.jl")
    include("utils.jl")

    function __init__()
        @require CUDA="052768ef-5323-5732-b1bb-66c8b64840ba" begin
            if CUDA.functional() && CUDA.has_cutensor()
                const CuArray = CUDA.CuArray
                const CuVector = CUDA.CuVector
                const CuMatrix = CUDA.CuMatrix
                const CuSVD = CUDA.CUSOLVER.CuSVD
                # scalar indexing is fine before 0.2
                # CUDA.allowscalar(false)
                include("cuda/base.jl") 
                include("cuda/contractions.jl")
                include("cuda/compressions.jl")
            end
        end
    end
end