"""
Functions that produce the "explicit" pretwitch + posttwitch dataframes
shared with collaborators — columns are
(lineage_name, minutes_post_first_cleavage, RL_micrometers,
VD_micrometers, AP_micrometers).

The RL axis is the flipped LR axis (RL = -LR), kept consistent with the
display/movie flip in `load_average_annotations`. The VD axis is **not**
negated for post-twitch data: the raw HDF5/`_for_ben` `z` is already in
the correct VD orientation. (Pre-twitch coords still need a DV flip; see
`get_pretwitch_explicit_df`.) Both phases emit the same `VD` column name
so the pre/post dataframes concatenate cleanly.

Pretwitch coordinates are taken from `ryan_data/final_20251205_pre-twitch_coords.csv`
via `get_pretwitch_explicit_df` (defined in `src/seam_cell_to_lineage_map.jl`).
Posttwitch coordinates come from a `_for_ben.csv` produced by
`resave_for_ben`, with the positional model cell names translated to
embryonic lineage names via the
`ryan_data/MIPAV_PositionalModel_Packer_Naming_Correlations_v6.csv` table.

`get_combined_explicit_df` returns the three dataframes; `write_combined_explicit_csvs`
also writes them to disk.
"""

using CSV: CSV
using DataFrames: DataFrame, ByRow, vcat, select, subset!
using Dates: Dates

"""
    get_seam_cells_explicit_df(avg_models; positional_to_lineage_dict)

Build the seam-cell rows for the posttwitch explicit dataframe. Reads
the seam cell points from each `avg_model` and translates their
positional names to lineage names via `positional_to_lineage_dict`.
Rows whose lineage couldn't be resolved are dropped at the caller.
"""
function get_seam_cells_explicit_df(
    avg_models;
    positional_to_lineage_dict::Dict,
)
    seam_cell_names = ["a0L", "a0R", "H0L", "H0R", "H1L", "H1R", "H2L", "H2R",
                       "V1L", "V1R", "V2L", "V2R", "V3L", "V3R", "V4L", "V4R",
                       "V5L", "V5R", "V6L", "V6R", "TL", "TR"]
    # Reorder right-then-left so the seam-cell-pts layout matches the
    # avg_model output.
    seam_cell_names = seam_cell_names[[2:2:end; 1:2:end]]
    seam_cell_lineage_names = get.(
        (positional_to_lineage_dict,),
        seam_cell_names,
        missing,
    )
    dfs = map(enumerate(avg_models)) do (i, model)
        pts = seam_cell_pts(model, 2)
        pts = swapyz_scale.(pts)
        # LR flipped: RL = -LR (-x). VD = +z (post-twitch z already in VD
        # orientation, not negated). AP (y) unchanged.
        DataFrame(
            lineage_name = seam_cell_lineage_names,
            minutes_post_first_cleavage = (i - 1) * 370 / (length(avg_models) - 1) + 381,
            RL_micrometers = pts .|> x -> -x[1],
            VD_micrometers = pts .|> x -> x[3],
            AP_micrometers = pts .|> x -> x[2],
        )
    end
    subset(vcat(dfs...), :lineage_name => ByRow(!ismissing))
end

"""
    get_combined_explicit_df(; avg_models, ben_csv_path,
                             pretwitch_df = get_pretwitch_df(),
                             annotation_name_translation_df = get_annotation_name_translation_df(),
                             add_unsmoothed_seam_cells = false)

Returns `(; pretwitch_explicit_df, posttwitch_for_ben_explicit_df, combined_df)`.

- `pretwitch_explicit_df`: pretwitch coords aligned to the first avg_model frame.
- `posttwitch_for_ben_explicit_df`: `_for_ben.csv` rows translated to
  the explicit schema, with positional names mapped to lineage names
  via `annotation_name_translation_df`.
- `combined_df`: vertical concatenation of the two.

The `_for_ben.csv` already contains the seam cells (added to the averaged
HDF5 as a `seam_cells` group and smoothed before export, then read back by
`resave_for_ben`). Setting `add_unsmoothed_seam_cells = true` *additionally*
appends the **unsmoothed** seam cells pulled directly from `avg_models` via
`get_seam_cells_explicit_df`. That double-counts every seam cell (a smoothed
copy and an unsmoothed copy per timepoint), so it defaults to `false`.
"""
function get_combined_explicit_df(;
    avg_models,
    ben_csv_path::AbstractString,
    pretwitch_df::DataFrame = get_pretwitch_df(),
    annotation_name_translation_df::DataFrame = get_annotation_name_translation_df(),
    add_unsmoothed_seam_cells::Bool = false,
)
    pretwitch_explicit_df = get_pretwitch_explicit_df(pretwitch_df; avg_models)

    positional_to_lineage_dict = Dict(
        annotation_name_translation_df.var"Positional Model Cell Name" .=>
        annotation_name_translation_df.var"Lineage Name",
    )

    ben_df = CSV.read(ben_csv_path, DataFrame)
    # LR flipped: RL = -LR (-x). VD = +z (post-twitch z already in VD
    # orientation, not negated). AP (y) unchanged.
    posttwitch_for_ben_explicit_df = select(ben_df,
        :cell => ByRow(cell -> get(positional_to_lineage_dict, cell, missing)) => :lineage_name,
        :time => :minutes_post_first_cleavage,
        :x => ByRow(-) => :RL_micrometers,
        :z => :VD_micrometers,
        :y => :AP_micrometers,
    )
    if add_unsmoothed_seam_cells
        posttwitch_for_ben_explicit_df = vcat(
            posttwitch_for_ben_explicit_df,
            get_seam_cells_explicit_df(avg_models; positional_to_lineage_dict),
        )
    end
    subset!(posttwitch_for_ben_explicit_df, :lineage_name => ByRow(!ismissing))

    combined_df = vcat(pretwitch_explicit_df, posttwitch_for_ben_explicit_df)

    return (;
        pretwitch_explicit_df,
        posttwitch_for_ben_explicit_df,
        combined_df,
    )
end

"""
    write_combined_explicit_csvs(; output_dir, avg_models, ben_csv_path,
                                   prefix = "", date_str = ...)

Write three CSVs into `output_dir` (with the default empty `prefix`):
  - `pretwitch_<date>.csv`
  - `posttwitch_<date>.csv`
  - `combined_<date>.csv`

A non-empty `prefix` is prepended as `<prefix>_…`.

Returns the same NamedTuple as `get_combined_explicit_df` with the
output paths attached.
"""
function write_combined_explicit_csvs(;
    output_dir::AbstractString,
    avg_models,
    ben_csv_path::AbstractString,
    prefix::AbstractString = "",
    date_str::AbstractString = Dates.format(Dates.now(), "yyyy_mm_dd_HHMMSS"),
    pretwitch_df::DataFrame = get_pretwitch_df(),
    annotation_name_translation_df::DataFrame = get_annotation_name_translation_df(),
    add_unsmoothed_seam_cells::Bool = false,
)
    mkpath(output_dir)
    result = get_combined_explicit_df(;
        avg_models = avg_models,
        ben_csv_path = ben_csv_path,
        pretwitch_df = pretwitch_df,
        annotation_name_translation_df = annotation_name_translation_df,
        add_unsmoothed_seam_cells = add_unsmoothed_seam_cells,
    )
    _p = isempty(prefix) ? "" : "$(prefix)_"
    pretwitch_path = joinpath(output_dir, "$(_p)pretwitch_$(date_str).csv")
    posttwitch_path = joinpath(output_dir, "$(_p)posttwitch_$(date_str).csv")
    combined_path = joinpath(output_dir, "$(_p)combined_$(date_str).csv")
    CSV.write(pretwitch_path, result.pretwitch_explicit_df)
    CSV.write(posttwitch_path, result.posttwitch_for_ben_explicit_df)
    CSV.write(combined_path, result.combined_df)
    return (; result..., pretwitch_path, posttwitch_path, combined_path)
end
