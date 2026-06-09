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
