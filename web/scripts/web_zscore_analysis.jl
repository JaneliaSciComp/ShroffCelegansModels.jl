# Thin wrapper around ShroffCelegansModelsWebInterface.ZscoreAnalysis.
# Data-loading (runtime includes) + server live in web/src/apps/ZscoreAnalysis.jl,
# whose precompile workload caches the Bonito.Table render path.
using ShroffCelegansModelsWebInterface: ZscoreAnalysis

if abspath(PROGRAM_FILE) == @__FILE__
    ZscoreAnalysis.main()
end
