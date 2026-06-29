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

include("apps/MeshscatterAverage.jl")

end # module ShroffCelegansModelsWebInterface
