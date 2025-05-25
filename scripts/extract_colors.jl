using CSV
using DataFrames

function extract_color(filename)
    df = CSV.read(filename, DataFrame; header=0)
    only(unique(eachrow(df[!,5:7])))
end
function extract_colors(
    dir = raw"\\nearline4.hhmi.org\shroff\shrofflab\model_building\For_Movies\0918_ColoredOutput"
)
    csv_files = filter(endswith(".csv"), readdir(dir))
    csv_files = joinpath.(dir, csv_files)
    colors = extract_color.(csv_files)
    stacked_colors = UInt8.(stack(colors; dims=1))
    stacked_colors_df = DataFrame(stacked_colors, ["R", "G", "B"])
    annotation_names = replace.(csv_files, ".csv" => "")
    name_df = DataFrame("Annotation" => basename.(annotation_names))
    return [name_df stacked_colors_df]
end
# CSV.write("colors.csv", extract_colors())