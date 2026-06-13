# Guards the flat step-6 priming refactor (prime_dataset_positions!). A full
# numerical-equivalence test needs the /nearline dataset+annotation data that only
# exists on the cluster (verified there via an HDF5 diff against a baseline). What we
# CAN check locally without that data: the function is callable with the documented
# signature, the empty-input path is a clean no-op, and the in-region BLAS thread
# count is restored afterward (the `finally` invariant — important because the pad
# loop pins BLAS to 1 while running).
using LinearAlgebra
using GeometryBasics: Point3

@testset "prime_dataset_positions! empty input + BLAS restore" begin
    lattice = joinpath(@__DIR__, "fixtures", "lattice.csv")
    model = ShroffCelegansModels.build_celegans_model(lattice)
    avg_models = [model, model, model]            # length must match timepoints
    timepoints = LinRange(0, 1, 3)

    cache = Dict{String, Vector{Vector{Point3{Float64}}}}()
    empty_datasets = Dict{String, Vector{ShroffCelegansModels.Datasets.NormalizedDataset}}()

    blas_before = BLAS.get_num_threads()
    grand_total = ShroffCelegansModels.prime_dataset_positions!(
        empty_datasets, cache, timepoints; avg_models,
    )
    @test grand_total == 0                         # no datasets ⇒ nothing primed
    @test isempty(cache)                           # cache untouched
    @test BLAS.get_num_threads() == blas_before    # restored in finally

    # Timepoints/avg_models length mismatch is rejected.
    @test_throws AssertionError ShroffCelegansModels.prime_dataset_positions!(
        empty_datasets, cache, LinRange(0, 1, 5); avg_models,
    )
end
