using HDF5
using Printf
using DataFrames
using CSV

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
            println("Processing strain: $strain")

            # Get annotation labels (cell names)
            labels = h5f[strain]["annotations"][]
            num_cells = length(labels)

            # Read all timepoints
            matrices = Matrix{Float64}[]
            for tp in 1:201
                tpk = @sprintf("timepoint_%03d", tp)
                ds = h5f[strain][tpk]
                push!(matrices, ds[])
            end

            # Stack into 3D array: [cells × timepoints × coordinates]
            positions = stack(matrices, dims=2)

            # Generate time vector (201 timepoints from start to end)
            start_time, end_time = time_range
            minutes = range(start_time, end_time, length=201) |> collect

            # Create rows for each cell × timepoint combination
            for (cell_idx, cell_name) in enumerate(labels)
                for (tp_idx, time) in enumerate(minutes)
                    # Extract x, y, z coordinates
                    x = positions[cell_idx, tp_idx, 1]
                    y = positions[cell_idx, tp_idx, 2]
                    z = positions[cell_idx, tp_idx, 3]

                    # Create row: cell, time, x, y, z
                    push!(all_rows, (cell=cell_name, time=time, x=x, y=y, z=z))
                end
            end
        end

        all_rows
    end

    println("Creating DataFrame with $(length(ben_data)) rows")

    # Convert to DataFrame
    df = DataFrame(ben_data)

    # Write to CSV (with index column as first column, unnamed)
    println("Writing CSV file: $target_filename")
    CSV.write(target_filename, df, writeheader=true)

    println("✓ Conversion complete!")
    println("  Rows: $(nrow(df))")
    println("  Unique cells: $(length(unique(df.cell)))")
    println("  Time range: $(extrema(df.time))")

    return df
end

# Example usage (commented out):
# resave_for_ben("edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5")
