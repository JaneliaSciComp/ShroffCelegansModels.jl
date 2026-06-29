# Thin wrapper around ShroffCelegansModelsWebInterface.ShowAverageAnnotations.
# Render/data logic does its PVC reads at runtime inside the submodule's main();
# shared render compilation is cached by the package (CommonScenes workload).
using ShroffCelegansModelsWebInterface: ShowAverageAnnotations

if abspath(PROGRAM_FILE) == @__FILE__
    ShowAverageAnnotations.main()
end
