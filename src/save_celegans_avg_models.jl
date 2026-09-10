using Dates
using Printf
using ShroffCelegansModels.HDF5

# needs modelio.jl
include("demo_averaging/modelio.jl")

avg_models_filename() = "celegans_avg_models_" * Dates.format(now(), "yyyy_mm_dd") * ".h5"

function save_avg_models(avg_models_filename = avg_models_filename(), avg_models = avg_models, offset = 0)
    h5open(avg_models_filename, "w") do h5f
        for (i, m) in pairs(avg_models)
            save_celegans_model(h5f, @sprintf("avg_model_%03d", i + offset), m)
        end
    end
end

function load_avg_models(avg_models_filename= avg_models_filename())
    avg_models = []
    h5open(avg_models_filename, "r") do h5f
        for k in keys(h5f)
            group_attrs = attrs(h5f[k])
            if haskey(group_attrs, "julia_type") &&
               contains(group_attrs["julia_type"], "ShroffCelegansModels.Types.CelegansModel")
                push!(avg_models, load_celegans_model(h5f[k]))
            end
        end
    end
    identity.(avg_models)
end

"""
    load_latest_avg_models(; dir=ENV["RECOMPUTE_OUTPUT_DIR"] or "/data/annotations/recompute", default_filename=nothing)

Load the recompute pipeline's average models — the newest `avg_models_n<N>.h5`
(by mtime) in the recompute output directory. The previously-bundled
`celegans_avg_models_*.h5` snapshot is outdated and no longer shipped, so there
is no built-in fallback; pass `default_filename` (and ensure it exists) for a
local/offline source.
"""
function load_latest_avg_models(;
    dir::AbstractString = get(ENV, "RECOMPUTE_OUTPUT_DIR", "/data/annotations/recompute"),
    default_filename::Union{Nothing,AbstractString} = nothing,
)
    if isdir(dir)
        candidates = String[
            joinpath(dir, f) for f in readdir(dir)
            if startswith(f, "avg_models_n") && endswith(f, ".h5") && !occursin(".tmp.", f)
        ]
        if !isempty(candidates)
            chosen = argmax(mtime, candidates)
            @info "Loading average models" chosen dir
            return load_avg_models(chosen)
        end
    end
    if default_filename !== nothing && isfile(default_filename)
        @info "Loading average models (fallback)" default_filename
        return load_avg_models(default_filename)
    end
    error("load_latest_avg_models: no avg_models_n*.h5 found in $dir and no usable fallback (default_filename=$default_filename)")
end

function save_measurements(avg_models_filename = avg_models_filename())
    h5open(avg_models_filename) do h5f
        attrs(h5f)["voxel_pitch_micrometers"] = 0.1625
        h5f["measurements/volume"] = ShroffCelegansModels.volume_by_cross_section.(avg_models) .* 0.1625^3
        get_length(model) = ShroffCelegansModels.central_spline(model)(1.0)[3]
        h5f["measurements/length"] = get_length.(avg_models) .* 0.1625
        attrs(h5f["measurements/volume"])["units"] = "micrometers^3"
        attrs(h5f["measurements/length"])["units"] = "micrometers"
    end
end