"""
Fetch the 72 per-dataset CSV archives from the OpenShift PVC and unpack them
into a dated folder on /nrs, ready for use with NEARLINE_BASE.

Target layout after this script:
  /nrs/shroff/data_internal/celegans_mipav_data/YYYY_MM_DD/shroff/shrofflab/...

Set NEARLINE_BASE=/nrs/shroff/data_internal/celegans_mipav_data/YYYY_MM_DD
before running the pipeline so all /nearline paths are remapped correctly.

Usage:
  julia scripts/unpack_csv_archives_to_nrs.jl [YYYY_MM_DD]

The date defaults to today. Run from the repo root.
"""

using Dates: today, format

const REPO_ROOT   = dirname(dirname(abspath(@__FILE__)))
const OC          = joinpath(REPO_ROOT, "deployment", ".pixi", "envs", "default", "bin", "oc")
const NAMESPACE   = "shroff-data-test"
const CONTAINER   = "show-average-annotations"
const ARCHIVE_DIR = "/data/annotations/csv_archives"
const NRS_BASE    = "/nrs/shroff/data_internal/celegans_mipav_data"

const ARCHIVES = [
    "CND-1_032423_Pos2.tar.gz",
    "CND-1_033023_Pos2.tar.gz",
    "CND-1_033023_Pos3.tar.gz",
    "DCR4221_040215_UTP 2.tar.gz",
    "DCR4221_040215_UTP 3.tar.gz",
    "DCR4221_040215_UTP 4.tar.gz",
    "DCR4221_040215_UTP 5.tar.gz",
    "DCR6485_RPM1_NU_011419_Pos0.tar.gz",
    "DCR6485_RPM1_NU_011419_Pos4.tar.gz",
    "DCR6485_RPM1_NU_021020_Pos2.tar.gz",
    "Efn-4_031723_Pos1.tar.gz",
    "Efn-4_031723_Pos3.tar.gz",
    "Efn-4_031723_Pos4.tar.gz",
    "JCC596_NU_082619_Pos3.tar.gz",
    "JCC596_NU_091119_Pos2.tar.gz",
    "JCC596_NU_091119_Pos3.tar.gz",
    "KP9305_NU_073019_Pos0.tar.gz",
    "KP9305_NU_073019_Pos2.tar.gz",
    "KP9305_NU_073019_Pos4.tar.gz",
    "OD1599_NU_112619_Pos0.tar.gz",
    "OD1599_NU_112719_Pos3.tar.gz",
    "OD1599_NU_120619_Pos2.tar.gz",
    "RW10131_052918_Pos0.tar.gz",
    "RW10131_052918_Pos1.tar.gz",
    "RW10131_retracked_202404_Pos1.tar.gz",
    "RW10131_retracked_202404_Pos4.tar.gz",
    "RW10131_retracked_202405_Pos1.tar.gz",
    "RW10375_Pos0.tar.gz",
    "RW10375_Pos1.tar.gz",
    "RW10375_Pos2.tar.gz",
    "RW10557_031021_Pos0.tar.gz",
    "RW10557_031021_Pos1.tar.gz",
    "RW10557_031021_Pos2.tar.gz",
    "RW10584_051817.tar.gz",
    "RW10584_052517.tar.gz",
    "RW10584_101017.tar.gz",
    "RW10598_Pos0.tar.gz",
    "RW10598_Pos1.tar.gz",
    "RW10598_Pos4.tar.gz",
    "RW10598_retracked_202307_Pos1.tar.gz",
    "RW10598_retracked_202307_Pos2.tar.gz",
    "RW10598_retracked_202307_Pos4.tar.gz",
    "RW10711_Pos0.tar.gz",
    "RW10711_Pos1.tar.gz",
    "RW10711_Pos3.tar.gz",
    "RW10742_Pos1.tar.gz",
    "RW10742_Pos4.tar.gz",
    "RW10742_Pos5.tar.gz",
    "RW10752_NU_022519_Pos0.tar.gz",
    "RW10752_NU_031219_Pos1.tar.gz",
    "RW10752_NU_031219_Pos2.tar.gz",
    "RW10753_Pos1.tar.gz",
    "RW10753_Pos2.tar.gz",
    "RW10753_Pos6.tar.gz",
    "RW10896_Pos1.tar.gz",
    "RW10896_Pos4.tar.gz",
    "RW10896_Pos6.tar.gz",
    "RW10896_retracked_202311_Pos1.tar.gz",
    "RW10896_retracked_202311_Pos2.tar.gz",
    "RW10896_retracked_202311_Pos3.tar.gz",
    "Vab-1_Pos0.tar.gz",
    "Vab-1_Pos2.tar.gz",
    "Vab-1_Pos3.tar.gz",
    "efn-1_Pos0.tar.gz",
    "efn-1_Pos3.tar.gz",
    "efn-1_Pos4.tar.gz",
    "efn-2_120222_Pos1.tar.gz",
    "efn-2_120222_Pos2.tar.gz",
    "efn-2_120222_Pos6.tar.gz",
    "lin-26_021323_Pos3.tar.gz",
    "lin-26_031523_Pos4.tar.gz",
    "lin-26_Pos4.tar.gz",
]

function find_running_pod()
    out = readchomp(`$OC get pod -n $NAMESPACE
        --field-selector=status.phase=Running
        -o jsonpath={.items[0].metadata.name}`)
    isempty(out) && error("No Running pod found in namespace $NAMESPACE")
    return out
end

function unpack_archive(pod::String, archive::String, target::String)
    # Stream: oc exec ... cat archive | tar xzf - -C target
    # Using pipeline so the compressed bytes never hit disk.
    cat_cmd = `$OC exec -n $NAMESPACE $pod -c $CONTAINER --
        sh -c "cat \"$ARCHIVE_DIR/$archive\""`
    tar_cmd = `tar xzf - -C $target`
    open(tar_cmd, "w") do tar_io
        open(cat_cmd) do cat_io
            while !eof(cat_io)
                write(tar_io, read(cat_io, 65536))
            end
        end
    end
end

function main()
    date_str = length(ARGS) >= 1 ? ARGS[1] : format(today(), "yyyy_mm_dd")
    target    = joinpath(NRS_BASE, date_str, "shroff")
    nearline_base = joinpath(NRS_BASE, date_str)

    pod = find_running_pod()
    println("Pod:          $pod")
    println("Target:       $target")
    println("NEARLINE_BASE will be: $nearline_base")
    println()

    mkpath(target)

    total = length(ARCHIVES)
    ok    = 0
    fail  = 0
    failed_archives = String[]

    for (i, archive) in enumerate(ARCHIVES)
        print("[$i/$total] $archive ... ")
        try
            unpack_archive(pod, archive, target)
            println("ok")
            ok += 1
        catch err
            println("FAILED: $err")
            fail += 1
            push!(failed_archives, archive)
        end
    end

    println()
    println("Done: $ok ok, $fail failed")
    if !isempty(failed_archives)
        println("Failed archives:")
        for a in failed_archives
            println("  $a")
        end
    end
    println()
    println("To use with the pipeline:")
    println("  export NEARLINE_BASE=\"$nearline_base\"")
    println("  julia --project=web web/scripts/run_recompute_if_needed.jl")
end

main()
