function get_seam_cells_explicit_df(avg_models=avg_models)
    seam_cell_names = ["a0L", "a0R", "H0L", "H0R", "H1L", "H1R", "H2L", "H2R", "V1L", "V1R", "V2L", "V2R", "V3L", "V3R", "V4L", "V4R", "V5L", "V5R", "V6L", "V6R", "TL", "TR"]
    # right then left
    seam_cell_names = seam_cell_names[[2:2:end; 1:2:end]]
    seam_cell_lineage_names = get.(
        (positional_to_lineage_dict,),
        seam_cell_names,
        missing
    )
    dfs = map(enumerate(avg_models)) do (i, model)
        pts = seam_cell_pts(model, 2)
        pts = swapyz_scale.(pts)
        DataFrame(
            lineage_name = seam_cell_lineage_names,
            minutes_post_first_cleavage = (i - 1) * 370 / (length(avg_models) - 1) +381, 
            LR_micrometers = pts .|> x -> x[1],
            DV_micrometers = pts .|> x -> x[3],
            AP_micrometers = pts .|> x -> x[2],
        )
    end
    subset(vcat(dfs...), :lineage_name => ByRow(!ismissing))
end 

pretwich_explicit_df = get_pretwitch_explicit_df()
CSV.write("2026_04_02_pretwitch.csv", pretwitch_explicit_df)

annotation_name_translation_df = get_annotation_name_translation_df()
positional_to_lineage_dict = Dict(annotation_name_translation_df.var"Positional Model Cell Name" .=> annotation_name_translation_df.var"Lineage Name")
df = CSV.read("smoothed_average_annotations_r020_theta020_z030_for_ben.csv", DataFrame)
posttwitch_for_ben_explicit_df = select(df, :cell => ByRow(cell -> positional_to_lineage_dict[cell]) => :lineage_name, :time => :minutes_post_first_cleavage, :x => :LR_micrometers, :z => :DV_micrometers, :y=> :AP_micrometers)
posttwitch_for_ben_explicit_df = vcat(posttwitch_for_ben_explicit_df, get_seam_cells_explicit_df(avg_models))
subset!(posttwitch_for_ben_explicit_df, :lineage_name => ByRow(!ismissing))
CSV.write("2026_04_02_posttwitch.csv", posttwitch_for_ben_explicit_df)




single_df = vcat(pretwitch_explicit_df, posttwitch_for_ben_explicit_df)

