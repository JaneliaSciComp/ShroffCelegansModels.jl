# Thin wrapper around ShroffCelegansModelsWebInterface.MeshscatterAverage2024,
# which reuses MeshscatterAverage's render/server logic (and precompile cache)
# with the 2024-10-24 data file, port 8590, and its own proxy path.
using ShroffCelegansModelsWebInterface: MeshscatterAverage2024

if abspath(PROGRAM_FILE) == @__FILE__
    MeshscatterAverage2024.main()
end
