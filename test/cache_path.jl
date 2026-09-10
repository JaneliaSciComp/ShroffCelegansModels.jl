# Tests for _latest_cache_path: it must prefer the recompute (PVC) file whenever
# it EXISTS, not by mtime — the image build resets the baked-in file's mtime to
# build time, so an mtime comparison would always pick the stale bundled copy.

@testset "_latest_cache_path prefers recompute file by existence, not mtime" begin
    fname = "cache_path_test_$(getpid()).h5"
    baked = joinpath(pkgdir(ShroffCelegansModels), fname)  # == joinpath(@__DIR__src, "..", fname)
    mktempdir() do recompute_dir
        recomputed = joinpath(recompute_dir, fname)
        try
            withenv("RECOMPUTE_OUTPUT_DIR" => recompute_dir) do
                # No recompute file -> fall back to the baked-in path.
                isfile(recomputed) && rm(recomputed)
                @test normpath(ShroffCelegansModels._latest_cache_path(fname)) == normpath(baked)

                # Recompute file present -> chosen even when baked-in is NEWER
                # (proves existence wins over mtime).
                touch(recomputed)
                sleep(0.05)
                touch(baked)   # baked-in now strictly newer than recomputed
                @test ShroffCelegansModels._latest_cache_path(fname) == recomputed
            end
        finally
            isfile(baked) && rm(baked)
        end
    end
end
