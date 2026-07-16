# Thin wrapper around ShroffCelegansModelsWebInterface.LatticeOrientation.
# The DOM/render logic and its PrecompileTools workload live in the package
# (web/src/apps/LatticeOrientation.jl).
using ShroffCelegansModelsWebInterface: LatticeOrientation

if abspath(PROGRAM_FILE) == @__FILE__
    LatticeOrientation.main()
end
