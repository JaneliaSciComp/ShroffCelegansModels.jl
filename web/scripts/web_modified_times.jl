# Thin wrapper around ShroffCelegansModelsWebInterface.ModifiedTimes.
# The DOM/render logic and its PrecompileTools workload live in the package
# (web/src/apps/ModifiedTimes.jl).
using ShroffCelegansModelsWebInterface: ModifiedTimes

if abspath(PROGRAM_FILE) == @__FILE__
    ModifiedTimes.main()
end
