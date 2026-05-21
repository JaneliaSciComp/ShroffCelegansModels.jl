"""
Scans every dataset listed in the active config_path and records the mtime of
each Decon_reg_*/.../integrated_annotation/annotations.csv file. Used by the
OpenShift CronJob to track when integrated annotations were last edited.

Writes two HDF5 files into the directory given by MODIFIED_TIMES_DIR
(default: /data/annotations/modified_times):

  - modified_times_YYYY_mm_dd_HHMMSS.h5  (timestamped snapshot, kept as history)
  - modified_times.h5                    (overwritten each run, "latest" pointer)
"""

using Dates: now, format
using ShroffCelegansModels
using ShroffCelegansModels.JSON3
using ShroffCelegansModels.MIPAVIO: get_modified_times_unix, save_modified_times_unix

include(joinpath(@__DIR__, "..", "..", "src", "demo_averaging", "read_config_json.jl"))

# Write `dst` atomically: invoke `f(tmp_path)` to produce the file, then `mv` it
# into place. POSIX rename() is atomic on the same filesystem, so readers always
# see either the previous version or the new one — never a partial write.
function write_atomic(f, dst::AbstractString)
    tmp = string(dst, ".tmp.", getpid(), ".", time_ns())
    try
        f(tmp)
        mv(tmp, dst; force=true)
    catch
        isfile(tmp) && rm(tmp; force=true)
        rethrow()
    end
end

function main()
    output_dir = get(ENV, "MODIFIED_TIMES_DIR", "/data/annotations/modified_times")
    mkpath(output_dir)

    config_path = ShroffCelegansModels.config_path
    @info "Loading datasets from config" config_path
    _, _, datasets = read_config_json(config_path)
    @info "Loaded datasets" groups=length(datasets) total=sum(length, values(datasets))

    @info "Scanning integrated_annotation mtimes"
    modified_times = get_modified_times_unix(datasets)

    timestamp = format(now(), "yyyy_mm_dd_HHMMSS")
    snapshot_path = joinpath(output_dir, "modified_times_$(timestamp).h5")
    latest_path = joinpath(output_dir, "modified_times.h5")

    @info "Writing snapshot" snapshot_path
    write_atomic(snapshot_path) do tmp
        save_modified_times_unix(datasets, modified_times; filepath=tmp)
    end

    @info "Updating latest" latest_path
    write_atomic(latest_path) do tmp
        save_modified_times_unix(datasets, modified_times; filepath=tmp)
    end

    @info "Done"
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
