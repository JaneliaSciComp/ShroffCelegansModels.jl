"""
    DebugApAxisRetrackLive

Submodule for the retracked live AP-axis debug viewer (deployment container
`debug-ap-axis-retrack-live`, port 9281). Like [`DebugApAxisLive`](@ref) but over
a hardcoded set of retracked datasets, with per-route error capture.

`NormalizedDataset` construction reads `cell_key` files from the PVC, so the
dataset list is built at **runtime** inside [`main`](@ref), not at module load.
"""
module DebugApAxisRetrackLive

using WGLMakie
using Bonito
using ShroffCelegansModels

module Launch
    using WGLMakie
    using Bonito
    using ShroffCelegansModels
end

const PORT = 9281
const PROXY_PATH = "debug_annotation_ap_axis_retrack_live"

_include_launch() = Base.include(Launch,
    joinpath(pkgdir(ShroffCelegansModels), "scripts", "launch_show_average_annotations.jl"))

# Built at runtime — the NormalizedDataset constructor reads cell_key files from disk.
function _retracked_datasets()
    d = Dict{String, Vector{ShroffCelegansModels.NormalizedDataset}}()
    d["RW10131_SLS268"] = ShroffCelegansModels.NormalizedDataset.([
        "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240401/RW10131_SLS6_New/Pos4/SPIMB/Reg_Sample/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240429/SLS268_RW10131_SLS6/Pos1/SPIMB/Reg_Sample/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10131/Data/2024_SLS268/20240507/Pos1/SPIMA/Reg_Sample/For_Tracking/RegB",
    ])
    d["RW10598"] = ShroffCelegansModels.NormalizedDataset.([
        "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos1/SPIMB/Reg_Sample/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230718/RW10598_NU/Pos2/SPIMB/Reg_Sample/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10598/2023_Data/Tracking/20230719/RW10598_NU/Pos4/SPIMB/Reg_Sample/For_Tracking/RegB",
    ])
    d["RW10896"] = ShroffCelegansModels.NormalizedDataset.([
        "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos1/SPIMB/Registered_Volumes/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos2/SPIMB/Registered_Volumes/For_Tracking/RegB",
        "/nearline/shroff/shrofflab/RW10896/Postwitching/2023_imaging/20231129/RW10896_NU/Pos3/SPIMB/Registered_Volumes/For_Tracking/RegB",
    ])
    return d
end

function _menu(datasets)
    DOM.div(DOM.ul(map(collect(keys(datasets))) do k
        DOM.li(k), DOM.ul(map(collect(keys(datasets[k]))) do i
            DOM.li(DOM.a(datasets[k][i].path, href="/$PROXY_PATH/$k/$i"))
        end)
    end))
end

function build_server(datasets)
    server = Server("0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/$PROXY_PATH/")
    route!(server, "/" => App(_menu(datasets)))
    for k in keys(datasets)
        for i in keys(datasets[k])
            route!(server, "/$k/$i" => App(; title="$k[$i]: Shroff C. elegans debug annotation AP axis RETRACK") do
                empty!(ShroffCelegansModels.annotations_cache)
                try
                    return Launch.debug_annotation_ap_axis(Launch.avg_models, datasets[k][i]; use_myuntwist=true)
                catch err
                    msg = sprint(showerror, err, catch_backtrace())
                    return DOM.pre(msg)
                end
            end)
        end
    end
    return server
end

function main()
    @info "Loading data (runtime include)"
    _include_launch()
    WGLMakie.activate!(; resize_to = :body)
    datasets = _retracked_datasets()
    @info "Launching server!" port=PORT
    server = build_server(datasets)
    if isinteractive()
        println("Press enter to quit")
        readline()
    else
        wait(server)
    end
end

end # module DebugApAxisRetrackLive
