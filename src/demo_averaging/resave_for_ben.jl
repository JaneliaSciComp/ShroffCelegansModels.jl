using HDF5
using Printf
using DataFrames
using CSV
using Statistics

"""
    resave_for_ben(filename; target_filename = replace(filename, ".h5" => "_for_ben.csv"), time_range=(381, 751))

Convert Mark's HDF5 format to Ben's CSV format for Transcriptome4D integration.

# Arguments
- `filename`: Input HDF5 file path (e.g., "edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5")
- `target_filename`: Output CSV file path (default: input filename with "_for_ben.csv" suffix)
- `time_range`: Tuple of (start_time, end_time) in minutes post first cleavage (mpfc)
  - Default: (381, 751) for post-twitch based on Ryan's new timeline
  - 381 mpfc: first post-twitch position
  - 751 mpfc: hatching time

# Output Format
CSV file with columns: cell, time, x, y, z
- One row per cell per timepoint
- Time in minutes post first cleavage (mpfc)

# New Timeline (Ryan's Email - Dec 2025)
- Dataset starts at 0 mpfc
- First positional data: 20 mpfc (four cell stage)
- Last pre-twitch: 380 mpfc
- First post-twitch: 381 mpfc
- Hatching: 751 mpfc

# Example
```julia
resave_for_ben("edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5")
# Output: edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells_for_ben.csv
# Time range: 381-751 mpfc (post-twitch)
```

# Notes
- Based on resave_for_tosif() function structure
- Uses Ryan's new timeline scale (mpfc instead of old mpf)
- Cell names from HDF5 annotation labels
- **Auto-detects number of timepoints** from HDF5 structure (no hard-coding)
- Works with any number of timepoints (201, 371, etc.)
"""
function resave_for_ben(filename;
                         target_filename = replace(filename, ".h5" => "_for_ben.csv"),
                         time_range = (381, 751))

    if isfile(target_filename)
        error("$target_filename exists")
    end

    println("Reading HDF5 file: $filename")

    # Extract data from HDF5
    ben_data = h5open(filename) do h5f
        all_rows = []

        for strain in keys(h5f)
            if strain == "timepoint_range"
                continue
            end
            println("Processing strain: $strain")

            # Get annotation labels (cell names)
            labels = h5f[strain]["annotations"][]
            num_cells = length(labels)

            # Auto-detect number of timepoints from HDF5 structure
            strain_keys = keys(h5f[strain])
            timepoint_keys = filter(k -> startswith(k, "timepoint_"), strain_keys)
            num_timepoints = length(timepoint_keys)
            println("  Detected $num_timepoints timepoints")

            # Read all timepoints
            matrices = Matrix{Float64}[]
            for tp in 1:num_timepoints
                tpk = @sprintf("timepoint_%03d", tp)
                ds = h5f[strain][tpk]
                push!(matrices, ds[])
            end

            # Stack into 3D array: [cells × timepoints × coordinates]
            positions = stack(matrices, dims=2)

            # Generate time vector (dynamically based on detected timepoints)
            start_time, end_time = time_range
            minutes = range(start_time, end_time, length=num_timepoints) |> collect
            println("  Time spacing: $(round((end_time - start_time) / (num_timepoints - 1), digits=4)) min/timepoint")

            # Create rows for each cell × timepoint combination
            for (cell_idx, cell_name) in enumerate(labels)
                for (tp_idx, time) in enumerate(minutes)
                    # Extract x, y, z coordinates
                    x = positions[cell_idx, tp_idx, 1]
                    y = positions[cell_idx, tp_idx, 2]
                    z = positions[cell_idx, tp_idx, 3]

                    # Create row: strain, cell, time, x, y, z
                    push!(all_rows, (strain=strain, cell=cell_name, time=time, x=x, y=y, z=z))
                end
            end
        end

        all_rows
    end

    println("Creating DataFrame with $(length(ben_data)) rows")

    # Convert to DataFrame (includes strain column)
    df_full = DataFrame(ben_data)

    # Find cells that appear in multiple strains
    cells_per_strain = combine(groupby(df_full, [:cell, :time]), :strain => (x -> length(unique(x))) => :n_strains)
    duplicate_cells = unique(cells_per_strain[cells_per_strain.n_strains .> 1, :cell])
    println("Found $(length(duplicate_cells)) cells appearing in multiple strains")

    # 1. For Ben: Averaged positions (cell, time, x, y, z)
    println("Creating averaged positions for Ben...")
    df_averaged = combine(groupby(df_full, [:cell, :time]),
                          :x => mean => :x,
                          :y => mean => :y,
                          :z => mean => :z)
    sort!(df_averaged, [:cell, :time])

    println("Writing averaged CSV for Ben: $target_filename")
    CSV.write(target_filename, df_averaged, writeheader=true)

    # 2. For Ryan: Raw data for duplicate cells only (strain, cell, time, x, y, z)
    ryan_duplicates_filename = replace(target_filename, ".csv" => "_ryan_duplicates.csv")
    df_duplicates = filter(row -> row.cell in duplicate_cells, df_full)
    sort!(df_duplicates, [:cell, :time, :strain])

    println("Writing duplicate cells CSV for Ryan: $ryan_duplicates_filename")
    CSV.write(ryan_duplicates_filename, df_duplicates, writeheader=true)

    # 3. For Ryan: Summary statistics (cell, time, mean_x, mean_y, mean_z, std_x, std_y, std_z, n_strains)
    ryan_stats_filename = replace(target_filename, ".csv" => "_ryan_stats.csv")
    df_stats = combine(groupby(df_full, [:cell, :time]),
                       :x => mean => :mean_x,
                       :y => mean => :mean_y,
                       :z => mean => :mean_z,
                       :x => std => :std_x,
                       :y => std => :std_y,
                       :z => std => :std_z,
                       :strain => (x -> length(unique(x))) => :n_strains)

    # Only keep rows where n_strains > 1
    df_stats = filter(row -> row.n_strains > 1, df_stats)
    sort!(df_stats, [:cell, :time])

    println("Writing statistics CSV for Ryan: $ryan_stats_filename")
    CSV.write(ryan_stats_filename, df_stats, writeheader=true)

    println("\n✓ Conversion complete!")
    println("  For Ben (averaged): $(nrow(df_averaged)) rows → $target_filename")
    println("  For Ryan (duplicates): $(nrow(df_duplicates)) rows → $ryan_duplicates_filename")
    println("  For Ryan (statistics): $(nrow(df_stats)) rows → $ryan_stats_filename")
    println("  Unique cells: $(length(unique(df_averaged.cell)))")
    println("  Cells in multiple strains: $(length(duplicate_cells))")
    println("  Time range: $(extrema(df_averaged.time))")

    return df_averaged
end

# Example usage (commented out):
# resave_for_ben("edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5")
