using Revise
using WGLMakie
using Bonito


if abspath(PROGRAM_FILE) == @__FILE__
    global run_web_main::Bool = true
end

include("../../scripts/launch_show_average_annotations.jl")

raw"""
RW10131_SLS268:
Y:\shrofflab\RW10131\Data\2024_SLS268\20240401\RW10131_SLS6_New\Pos4\SPIMB\Reg_Sample\For_Tracking
Y:\shrofflab\RW10131\Data\2024_SLS268\20240429\SLS268_RW10131_SLS6\Pos1\SPIMB\Reg_Sample\For_Tracking
Y:\shrofflab\RW10131\Data\2024_SLS268\20240507\Pos1\SPIMA\Reg_Sample\For_Tracking

RW10598:
Y:\shrofflab\RW10598\2023_Data\Tracking\20230718\RW10598_NU\Pos1\SPIMB\Reg_Sample\For_Tracking
Y:\shrofflab\RW10598\2023_Data\Tracking\20230718\RW10598_NU\Pos2\SPIMB\Reg_Sample\For_Tracking
Y:\shrofflab\RW10598\2023_Data\Tracking\20230719\RW10598_NU\Pos4\SPIMB\Reg_Sample\For_Tracking

RW10896:
Y:\shrofflab\RW10896\Postwitching\2023_imaging\20231129\RW10896_NU\Pos1\SPIMB\Registered_Volumes\For_Tracking
Y:\shrofflab\RW10896\Postwitching\2023_imaging\20231129\RW10896_NU\Pos2\SPIMB\Registered_Volumes\For_Tracking
Y:\shrofflab\RW10896\Postwitching\2023_imaging\20231129\RW10896_NU\Pos3\SPIMB\Registered_Volumes\For_Tracking
"""

"""
    datasets = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
    map(config_json.data.strains) do strain
        datasets[strain.name] = map(strain.folderpaths) do folder_path
            ShroffCelegansModels.NormalizedDataset(joinpath(folder_path, "RegB"))
        end
    end
"""

const retracked_datasets = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
retracked_datasets["RW10131_SLS268"] = ShroffCelegansModels.NormalizedDataset.([
    "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240401/RW10131_SLS6_New/Pos4/SPIMB/Reg_Sample/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240429/SLS268_RW10131_SLS6/Pos1/SPIMB/Reg_Sample/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240507/Pos1/SPIMA/Reg_Sample/For_Tracking/RegB"
])
retracked_datasets["RW10598"] = ShroffCelegansModels.NormalizedDataset.([
    "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos1/SPIMB/Reg_Sample/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos2/SPIMB/Reg_Sample/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230719/RW10598_NU/Pos4/SPIMB/Reg_Sample/For_Tracking/RegB"
])
retracked_datasets["RW10896"] = ShroffCelegansModels.NormalizedDataset.([
    "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos1/SPIMB/Registered_Volumes/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos2/SPIMB/Registered_Volumes/For_Tracking/RegB",
    "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos3/SPIMB/Registered_Volumes/For_Tracking/RegB"
])

function web_debug_annotation_ap_axis(datasets = retracked_datasets)
    menu = DOM.div(
        DOM.ul(
            map(collect(keys(datasets))) do k
                DOM.li(k),
                DOM.ul(
                    map(collect(keys(datasets[k]))) do i
                        DOM.li(DOM.a(datasets[k][i].path, href="/debug_annotation_ap_axis_retrack_live/$k/$i"))
                    end
                )
            end
        )
    )
    server = Server(
        "shroff-data.int.janelia.org", 9281;
        proxy_url="https://shroff-data.int.janelia.org/debug_annotation_ap_axis_retrack_live/"
    )
    route!(server, "/" => App(menu))
    for k in keys(datasets)
        for i in keys(datasets[k])
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans debug annotation AP axis RETRACK") do
                empty!(annotations_cache)
                #empty!(annotation_position_cache)
                #empty!(my_annotation_position_cache)
                try
                    return debug_annotation_ap_axis(avg_models, datasets[k][i]; use_myuntwist=true);
                catch err
                    msg = sprint(showerror, err, catch_backtrace())
                    return DOM.pre(msg)
                end
            end)
        end
    end
    return server
end

function web_main()
    WGLMakie.activate!(; resize_to = :body)
    server = web_debug_annotation_ap_axis()
end

println(@__FILE__)
println(abspath(PROGRAM_FILE))

if run_web_main
    wait(web_main())
    println("Press enter to quit")
    readline()
end
