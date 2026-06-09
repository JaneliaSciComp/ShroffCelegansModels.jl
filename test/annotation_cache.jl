# Tests for the annotation-cache loaders:
#  - load_annotation_cache reconstructs keys root-aware (Windows drive letter vs
#    Linux "nearline" root written by the recompute pipeline), returns the cache,
#    and is graceful when the file is missing.
#  - load_annotations_cache reads per-timepoint "mtime" attrs (the read(::Float64)
#    bug) and is graceful when missing.
#  - load_average_annotations is graceful when missing.

using ShroffCelegansModels.HDF5

@testset "load_annotation_cache root-aware keys + return + graceful" begin
    C = ShroffCelegansModels.my_annotation_position_cache
    mktempdir() do dir
        f = joinpath(dir, "mypos.h5")
        h5open(f, "w") do h5
            # A Windows-drive root (single-char "X") and a Linux "nearline" root.
            # Datasets are NxN matrices; eachrow -> Point3.
            h5["X/foo/bar/001"] = rand(2, 3)
            h5["X/foo/bar/002"] = rand(2, 3)
            h5["nearline/shroff/baz/001"] = rand(2, 3)
            h5["nearline/shroff/baz/002"] = rand(2, 3)
        end
        empty!(C)
        ret = ShroffCelegansModels.load_annotation_cache(filename = f)
        @test ret === C                                      # returns the populated cache
        # Linux "nearline" root reconstructed as an absolute path (== dataset.path).
        @test haskey(C, "/nearline/shroff/baz")
        @test length(C["/nearline/shroff/baz"]) == 2         # two timepoints
        # Single-char root kept in Windows "X:\…" form (later mapped by alias).
        @test any(k -> startswith(k, "X:") && occursin("foo", k), keys(C))
    end

    empty!(C)
    missing_file = joinpath(tempdir(), "missing_mypos_$(getpid()).h5")
    @test ShroffCelegansModels.load_annotation_cache(filename = missing_file) === C
    @test isempty(C)
end

@testset "load_annotations_cache reads finite mtime attrs + graceful" begin
    AC = ShroffCelegansModels.annotations_cache
    mktempdir() do dir
        f = joinpath(dir, "annots.h5")
        h5open(f, "w") do h5
            g = create_group(h5, "X/foo/bar")
            attrs(g)["range_start"] = 1
            attrs(g)["range_end"] = 2
            attrs(g)["new_untwist"] = UInt8(1)
            for (idx, mt) in ((1, 111.5), (2, 222.5))
                gd = create_group(g, lpad(idx, 3, '0'))
                attrs(gd)["mtime"] = mt          # Float64 attr (this is what broke read())
                write_dataset(gd, "cellA", [1.0, 2.0, 3.0])
            end
        end
        empty!(AC)
        ShroffCelegansModels.load_annotations_cache(filename = f)
        vals = collect(values(AC))
        @test length(vals) == 1
        # dataset-level mtime = max of per-idx finite mtimes (read without error).
        @test vals[1].mtime == 222.5
    end

    empty!(AC)
    missing_file = joinpath(tempdir(), "missing_annots_$(getpid()).h5")
    @test ShroffCelegansModels.load_annotations_cache(filename = missing_file) === AC
    @test isempty(AC)
end

@testset "load_average_annotations graceful on missing file" begin
    d = ShroffCelegansModels.load_average_annotations(
        filename = joinpath(tempdir(), "missing_avg_$(getpid()).h5"),
    )
    @test d isa Dict
    @test isempty(d)
end
