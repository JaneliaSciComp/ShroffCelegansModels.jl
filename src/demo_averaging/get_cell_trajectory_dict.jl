using BSplineKit

function get_cell_trajectory_dict(dataset::ShroffCelegansModels.Datasets.NormalizedDataset; use_myuntwist = false)
    annotations = load_straightened_annotations_over_time(dataset; use_myuntwist)
    cells = collect(keys(annotations[1]))

    cell_trajectory_dict = map(cells) do cell
        trajectory = map(eachindex(annotations)) do i
            if ismissing(annotations[i]) || !haskey(annotations[i], cell) || any(!isfinite, annotations[i][cell])
                missing
            else
                annotations[i][cell]
            end
        end
        pairs = Iterators.filter(((i,p),)->!ismissing(p),zip(eachindex(trajectory), trajectory))
        pairs = collect(pairs)
        # cell => (first.(pairs), last.(pairs))
        try
            if length(pairs) < 2
                cell => missing
            else
                cell => extrapolate(
                    interpolate((first.(pairs) .- 1)/(length(trajectory)-1), last.(pairs), BSplineKit.BSplineOrder(2)),
                    Flat()
                )
            end
        catch err
            @error "There was an issue getting cell trajectory" cell dataset.path trajectory
            rethrow()
        end
    end |> Dict

    return cell_trajectory_dict
end