using ShroffCelegansModels
using ShroffCelegansModels: get_pretwitch_df, get_pretwitch_points_at_time,
    pretwitch_reference_axes, pretwitch_axes_at_time, pretwitch_cells_by_name,
    pretwitch_lineage_track, Point3f, normalize, dot
using Statistics: mean, std

# Approach 3 (pretwitch_orientation.txt item 3): treat Approach 1's WormGuides
# compass interpolation as a working hypothesis ("expert knowledge") and look
# for individual pretwitch cell lineages whose position relative to the body
# centroid stays consistently aligned with Approach 1's rotating DV/LR axes
# across the whole pretwitch window. Also specifically "work back" (via full
# ancestor-lineage reconstruction) from cells already identified as DV/LR
# markers in the posttwitch lattice_orientation.jl / check_lattice_orientation.jl
# QC (Cpaaaa, ABplaappap, Caaaaa for DV; ABarpaappp for LR) to see whether
# their pretwitch ancestry also tracks Approach 1's axes.

pretwitch_df = get_pretwitch_df()
ref = pretwitch_reference_axes(pretwitch_df)
grouped = pretwitch_cells_by_name(pretwitch_df)

times = sort(unique(pretwitch_df.time))

println("Precomputing per-frame centroid + Approach-1 axes for $(length(times)) timepoints...")
centroid_of = Dict{Int,Point3f}()
axes_of = Dict{Int,Any}()
for t in times
    pts = collect(values(get_pretwitch_points_at_time(pretwitch_df, t)))
    centroid_of[t] = Point3f(sum(pts) / length(pts))
    axes_of[t] = pretwitch_axes_at_time(ref, t)
end

const MIN_DIST = 10f0  # ignore frames where the cell is too close to centroid for direction to be meaningful

"""
    axis_alignment(name)

Reconstruct `name`'s full ancestor-lineage track and return
`(n_frames, mean_dv, std_dv, mean_lr, std_lr)` scoring how consistently its
direction from the per-frame body centroid aligns with Approach 1's `dv(t)`
and `lr(t)` axes.
"""
function axis_alignment(name)
    ts, pts = pretwitch_lineage_track(name, grouped)
    dv_scores = Float64[]
    lr_scores = Float64[]
    for (t, p) in zip(ts, pts)
        v = p - centroid_of[t]
        norm(v) < MIN_DIST && continue
        vhat = normalize(v)
        a = axes_of[t]
        push!(dv_scores, dot(vhat, a.dv))
        push!(lr_scores, dot(vhat, a.lr))
    end
    isempty(dv_scores) && return (0, NaN, NaN, NaN, NaN)
    return (length(dv_scores), mean(dv_scores), std(dv_scores), mean(lr_scores), std(lr_scores))
end

norm(v) = sqrt(sum(abs2, v))

println()
println("=== Cross-check: known posttwitch DV/LR marker cells (lattice_orientation.jl), tracked back through pretwitch ===")
known = [
    ("Cpaaaa", :dv, "hyp7_Cpaaaa (the original DV-swap check cell)"),
    ("ABplaappap", :dv, "hyp6_ABplaappap (DV survey candidate)"),
    ("Caaaaa", :dv, "hyp7_Caaaaa (DV survey candidate)"),
    ("ABarpaappp", :lr, "hyp7_ABarpaappp (LR survey candidate)"),
]
for (name, axis, label) in known
    n, mdv, sdv, mlr, slr = axis_alignment(name)
    println("$label:")
    println("  n_frames=$n  mean_dv=$(round(mdv,digits=3)) std_dv=$(round(sdv,digits=3))  mean_lr=$(round(mlr,digits=3)) std_lr=$(round(slr,digits=3))")
end

println()
println("=== Broad scan: every leaf lineage present at the final pretwitch timepoint ===")
last_t = maximum(times)
leaves = collect(keys(get_pretwitch_points_at_time(pretwitch_df, last_t)))
println("Scanning $(length(leaves)) leaf lineages...")

results = map(leaves) do name
    n, mdv, sdv, mlr, slr = axis_alignment(name)
    (name=name, n=n, mean_dv=mdv, std_dv=sdv, mean_lr=mlr, std_lr=slr)
end
results = filter(r -> r.n >= 100, results)  # require a long, statistically meaningful track

println("($(length(results)) of $(length(leaves)) leaves have >=100 aligned frames)")

println()
println("--- Top 15 DV-axis candidates (by |mean_dv|, low std_dv preferred) ---")
dv_ranked = sort(results, by=r -> -abs(r.mean_dv) + r.std_dv)
for r in first(dv_ranked, 15)
    println("$(r.name): n=$(r.n) mean_dv=$(round(r.mean_dv,digits=3)) std_dv=$(round(r.std_dv,digits=3))")
end

println()
println("--- Top 15 LR-axis candidates (by |mean_lr|, low std_lr preferred) ---")
lr_ranked = sort(results, by=r -> -abs(r.mean_lr) + r.std_lr)
for r in first(lr_ranked, 15)
    println("$(r.name): n=$(r.n) mean_lr=$(round(r.mean_lr,digits=3)) std_lr=$(round(r.std_lr,digits=3))")
end
