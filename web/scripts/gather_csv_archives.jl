"""
For each dataset listed in the active config, walks the dataset's
`/nearline` path and tar+gzips every `.csv` file (preserving the
directory structure relative to `/nearline/shroff/`) into
`<CSV_ARCHIVE_DIR>/<group>_<date_position_tag>.tar.gz`.

The output dir defaults to `/data/annotations/csv_archives` so the
archives land on the same PVC as `/data/annotations/recompute/` and
can be served by nginx (`location /csv_archives/`).

Idempotent: skips a dataset if its archive file already exists; pass
`OVERWRITE=1` to force re-creation.

Invoke via a K8s Job (see deployment/shroff-data-test/10-job-csv-archives.yaml).
"""

using ShroffCelegansModels
using ShroffCelegansModels: read_config_json, date_position_tag

const NEARLINE_BASE = "/nearline/shroff"

"""
    to_unix_nearline(path)

Convert a dataset.path that may be Windows-style (`X:\\shrofflab\\...`)
into the Linux-mounted `/nearline/shroff/...` form. Mirrors the
conversion in src/recompute_pipeline.jl and src/parse_worm_dataset_path.jl.
"""
function to_unix_nearline(p::AbstractString)
    s = String(p)
    if Sys.isunix() && occursin('\\', s)
        if length(s) >= 3 && isuppercase(s[1]) && s[2] == ':' && s[3] == '\\'
            s = NEARLINE_BASE * "/" * s[4:end]
        end
        s = replace(s, "\\" => "/")
    end
    return s
end

"""
    gather_dataset_csvs(dataset_path_linux, archive_path)

Run `find <rel> -name "*.csv" -print0 | tar czf <archive> -C /nearline/shroff
--null --files-from -`, so archive contents have paths relative to
`/nearline/shroff/`. Unpacking into a directory tree restores the
shroff-relative structure.
"""
function gather_dataset_csvs(
    dataset_path_linux::AbstractString,
    archive_path::AbstractString;
    base::AbstractString = NEARLINE_BASE,
)
    isdir(dataset_path_linux) || error("Dataset path does not exist: $dataset_path_linux")
    rel = relpath(dataset_path_linux, base)
    tmp = string(archive_path, ".tmp.", getpid(), ".", time_ns())
    try
        # Cmd(..., dir=base) makes `find rel ...` resolve relative to base
        # without changing process cwd.
        find_cmd = Cmd(`find $rel -type f -name "*.csv" -print0`; dir = base)
        tar_cmd  = `tar czf $tmp -C $base --null --files-from -`
        run(pipeline(find_cmd, tar_cmd))
        mv(tmp, archive_path; force = true)
    catch
        isfile(tmp) && rm(tmp; force = true)
        rethrow()
    end
    return archive_path
end

function main()
    output_dir = get(ENV, "CSV_ARCHIVE_DIR", "/data/annotations/csv_archives")
    overwrite = get(ENV, "OVERWRITE", "0") == "1"
    mkpath(output_dir)

    @info "Gathering CSV archives" output_dir overwrite
    _, _, datasets = read_config_json()
    flattened = collect((group, ds) for (group, dsets) in datasets for ds in dsets)
    n_total = length(flattened)
    n_done = 0
    n_skipped = 0
    n_failed = 0
    started_at = time()

    for (i, (group, ds)) in enumerate(flattened)
        tag = try
            date_position_tag(ds.path)
        catch
            "unparsed-$(i)"
        end
        archive_id = string(group, "_", tag)
        archive_path = joinpath(output_dir, "$(archive_id).tar.gz")

        if isfile(archive_path) && !overwrite
            n_skipped += 1
            @info "Skipping (exists)" archive_id archive_path
            continue
        end

        ds_path_linux = to_unix_nearline(ds.path)
        t0 = time()
        try
            gather_dataset_csvs(ds_path_linux, archive_path)
            sz = filesize(archive_path)
            n_done += 1
            @info "Archived" idx=i total=n_total archive_id size_bytes=sz elapsed_s=round(time() - t0; digits=2) path=archive_path
        catch err
            n_failed += 1
            @warn "Archive failed" archive_id ds_path_linux err
        end
    end

    @info "CSV archive gather complete" n_done n_skipped n_failed n_total elapsed_s=round(time() - started_at; digits=2)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
