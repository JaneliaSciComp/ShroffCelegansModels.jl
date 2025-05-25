using GLMakie
include("../../scripts/launch_show_average_annotations.jl")

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end