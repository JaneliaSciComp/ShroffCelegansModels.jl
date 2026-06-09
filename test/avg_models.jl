# Tests for load_latest_avg_models: it loads the newest avg_models_n<N>.h5 from
# the recompute output dir (RECOMPUTE_OUTPUT_DIR), and errors when none exist and
# no fallback is given. Round-trips a real model built from the lattice fixture.

@testset "load_latest_avg_models picks newest avg_models_n*.h5" begin
    lattice = joinpath(@__DIR__, "fixtures", "lattice.csv")
    model = ShroffCelegansModels.build_celegans_model(lattice)

    mktempdir() do dir
        # Older file with 1 model, then a newer file with 2 models.
        ShroffCelegansModels.save_avg_models(joinpath(dir, "avg_models_n5.h5"), [model])
        sleep(0.05)
        ShroffCelegansModels.save_avg_models(joinpath(dir, "avg_models_n7.h5"), [model, model])
        # A non-matching .h5 should be ignored.
        ShroffCelegansModels.save_avg_models(joinpath(dir, "something_else.h5"), [model])

        withenv("RECOMPUTE_OUTPUT_DIR" => dir) do
            m = ShroffCelegansModels.load_latest_avg_models()
            @test length(m) == 2                       # chose the newest (2-model) file
            @test all(x -> x isa ShroffCelegansModels.CelegansModel, m)
        end
    end

    # No avg_models_n*.h5 anywhere and no fallback -> error.
    mktempdir() do empty_dir
        withenv("RECOMPUTE_OUTPUT_DIR" => empty_dir) do
            @test_throws Exception ShroffCelegansModels.load_latest_avg_models()
        end
    end
end
