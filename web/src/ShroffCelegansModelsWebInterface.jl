"""
    ShroffCelegansModelsWebInterface

Package that hosts the Shroff Lab C. elegans Bonito/WGLMakie web applications.

Each web application lives in its **own submodule** (under `src/apps/`) for
isolation. Every submodule carries a [`PrecompileTools`](@ref) workload that
exercises its render and Bonito static-HTML serialize path against tiny
synthetic data, so the expensive Makie → WGLMakie → Bonito compilation is baked
into this package's precompile image instead of being repeated on every
container cold start.

Each `web/scripts/web_*.jl` entry point is now a thin wrapper that calls the
corresponding submodule's `main()`.
"""
module ShroffCelegansModelsWebInterface

# Shared render-primitive workload (helps the launch-script apps, whose own code
# can't run at build time).
include("CommonScenes.jl")

# One submodule per web app. MeshscatterAverage must come before its 2024 variant.
include("apps/MeshscatterAverage.jl")
include("apps/MeshscatterAverage2024.jl")
include("apps/ModifiedTimes.jl")
include("apps/LatticeOrientation.jl")
include("apps/ZscoreAnalysis.jl")
include("apps/ShowAverageAnnotations.jl")
include("apps/DebugApAxis.jl")
include("apps/DebugApAxisLive.jl")
include("apps/DebugApAxisRetrackLive.jl")
include("apps/FixApAxis.jl")

end # module ShroffCelegansModelsWebInterface
