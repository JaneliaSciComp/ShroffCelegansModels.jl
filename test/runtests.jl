using ShroffCelegansModels
using Test

@testset "ShroffCelegansModels.jl" begin
    include("cache_path.jl")
    include("annotation_cache.jl")
    include("avg_models.jl")
end
