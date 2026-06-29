# Thin wrapper around ShroffCelegansModelsWebInterface.DebugApAxisLive.
using ShroffCelegansModelsWebInterface: DebugApAxisLive

if abspath(PROGRAM_FILE) == @__FILE__
    DebugApAxisLive.main()
end
