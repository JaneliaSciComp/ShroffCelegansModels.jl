# Tests for the video helpers in scripts/crop_video.jl. These need VideoIO /
# ColorTypes / Makie, which live in the web environment, so run them with:
#
#     julia --project=web web/test/runtests.jl
#
using Test
using ColorTypes
const N0f8 = ColorTypes.N0f8

# Loads VideoIO, ColorTypes, Makie (all in the web env) and defines
# sane_framerate, crop_bounds, crop_video.
include(joinpath(@__DIR__, "..", "..", "scripts", "crop_video.jl"))

# Defines annotation_changes_pending (main() is guarded, so it doesn't run).
include(joinpath(@__DIR__, "..", "scripts", "run_recompute_if_needed.jl"))

@testset "crop_video helpers" begin
    @testset "sane_framerate normalizes bad rates" begin
        # Makie-saved VideoStreams report 1//0, which libx264 rejects (EINVAL -22).
        @test sane_framerate(1 // 0) == 24
        @test sane_framerate(0) == 24
        @test sane_framerate(-5) == 24
        @test sane_framerate(NaN) == 24
        # Valid rates pass through unchanged.
        @test sane_framerate(24) == 24
        @test sane_framerate(30 // 1) == 30 // 1
        @test sane_framerate(1 // 0; default = 15) == 15
    end

    @testset "crop_bounds yields even dimensions" begin
        img = fill(RGB{N0f8}(0, 0, 0), 200, 200)
        # odd-length white region (71 x 65); crop_bounds must round up to even
        # (yuv420p / libx264 require even width and height).
        img[60:130, 70:134] .= RGB{N0f8}(1, 1, 1)
        b = crop_bounds(img; offset = 50)
        @test all(r -> iseven(length(r)), b)
    end
end

@testset "annotation_changes_pending" begin
    mktempdir() do d
        changes = joinpath(d, "annotation_changes.h5")
        outdir = joinpath(d, "recompute"); mkpath(outdir)
        withenv("ANNOTATION_CHANGES_PATH" => changes, "RECOMPUTE_OUTPUT_DIR" => outdir) do
            @test annotation_changes_pending() == false          # no edits file yet
            out = joinpath(outdir, "avg_models_n5.h5")
            touch(out); sleep(0.05); touch(changes); sleep(0.05); touch(out)
            @test annotation_changes_pending() == false          # last run newer than edits
            sleep(0.05); touch(changes)
            @test annotation_changes_pending() == true           # edits newer than last run
        end
    end
    # No prior recompute output at all -> pending.
    mktempdir() do d
        changes = joinpath(d, "annotation_changes.h5"); touch(changes)
        withenv("ANNOTATION_CHANGES_PATH" => changes,
                "RECOMPUTE_OUTPUT_DIR" => joinpath(d, "empty")) do
            @test annotation_changes_pending() == true
        end
    end
end
