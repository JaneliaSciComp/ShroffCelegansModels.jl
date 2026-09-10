# Thin wrapper around ShroffCelegansModelsWebInterface.MeshscatterAverage.
# All render/server logic and the PrecompileTools workload live in the package
# (web/src/apps/MeshscatterAverage.jl) so the expensive Makie → WGLMakie → Bonito
# compilation is baked into the precompile image and skipped on cold start.
using ShroffCelegansModelsWebInterface: MeshscatterAverage

if abspath(PROGRAM_FILE) == @__FILE__
    MeshscatterAverage.main()
end
