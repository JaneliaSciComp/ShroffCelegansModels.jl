using ThinPlateSplines: tps_solve, tps_solve!

# `ws` (optional): a `Ref` holding a reusable `ThinPlateSplines.TPSWorkspace`, used
# to avoid re-allocating the ~tens-of-MB TPS solve buffers on every call when
# warping repeatedly against a fixed-size control lattice (the recompute averaging
# loop solves ~26k times against the constant 41×32 grid). The workspace is
# allocated lazily on first use (so its size matches the actual lattice) and reused
# thereafter. `ws === nothing` (the default) preserves the original allocating path.
# A `Ref` must not be shared across threads — allocate one per task.
function transform_annotations(from_model, to_model, annotations::AbstractVector{<: Point}; ws=nothing)
    warp_from = ShroffCelegansModels.lattice(from_model) |> vec
    warp_to = ShroffCelegansModels.lattice(to_model) |> vec
    tps_solved = if ws === nothing
        tps_solve(warp_from, warp_to, 1)
    else
        if ws[] === nothing
            P = eltype(warp_from)
            ws[] = ThinPlateSplines.TPSWorkspace{eltype(P)}(length(warp_from), length(P))
        end
        tps_solve!(ws[], warp_from, warp_to, 1)
    end

    transformed_pts = ThinPlateSplines.tps_deform(annotations, tps_solved)
    return transformed_pts
end
function transform_annotations(from_model, to_model, annotations::Dict; ws=nothing)
    return Dict(keys(annotations) .=> transform_annotations(from_model, to_model, collect(values(annotations)); ws))
end