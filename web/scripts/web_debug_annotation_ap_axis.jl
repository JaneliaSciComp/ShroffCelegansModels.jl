# Thin wrapper around ShroffCelegansModelsWebInterface.DebugApAxis.
using ShroffCelegansModelsWebInterface: DebugApAxis

if abspath(PROGRAM_FILE) == @__FILE__
    DebugApAxis.main()
end
