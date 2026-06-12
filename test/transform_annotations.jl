# Regression test for the annotation-averaging TPS warp.
#
# `transform_annotations` passes `lattice(model) |> vec` (a Vector{Point3}) to
# ThinPlateSplines' `tps_solve`/`tps_deform`. The core fork methods require a
# numeric K×D matrix and throw `BoundsError: ... Tuple{Int64} at index [2]` on a
# 1-D Point vector — which the averaging step silently swallowed, filling every
# non-seam annotation with Point3(NaN) (the recompute movie showed only seam
# cells). The `mkitti-geometrybasics-ext` branch adds a package extension whose
# `tps_solve`/`tps_deform` accept `AbstractVector{<:Point}`. Pin to it and this
# whole path works on Point vectors again.

using GeometryBasics: Point
import ThinPlateSplines

@testset "transform_annotations warps Point vectors (GeometryBasics ext)" begin
    # The GeometryBasics extension must actually be loaded, else tps_solve only
    # accepts numeric matrices and transform_annotations throws on Point input.
    ext = Base.get_extension(ThinPlateSplines, :ThinPlateSplinesGeometryBasicsExt)
    @test ext !== nothing

    lattice_csv = joinpath(@__DIR__, "fixtures", "lattice.csv")
    model = ShroffCelegansModels.build_celegans_model(lattice_csv)

    # A few real annotation positions taken from the model's own lattice.
    annotations = vec(ShroffCelegansModels.lattice(model))[1:5]
    @test eltype(annotations) <: Point

    # Identity warp (from == to): the result must be finite and reproduce the
    # input points. Before the fix this threw a BoundsError instead.
    warped = ShroffCelegansModels.transform_annotations(model, model, annotations)

    @test length(warped) == length(annotations)
    @test all(p -> all(isfinite, p), warped)
    @test all(i -> isapprox(warped[i], annotations[i]; atol = 1e-6), eachindex(annotations))
end
