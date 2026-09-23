using ShroffCelegansModels
using ShroffCelegansModels: get_pretwitch_df, get_pretwitch_points_at_time,
    pretwitch_reference_axes, pretwitch_axes_at_time, pretwitch_cells_by_name,
    Point3f, normalize, dot
using LinearAlgebra: I as LinearAlgebraI
using Statistics: mean

# v2: the first attempt (regress_pretwitch_axis_weights.jl) used each LEAF
# lineage's full ancestor-spliced track as a feature. That's badly collinear:
# sibling leaves that only diverged in the last few frames share an IDENTICAL
# trajectory for most of development, so ridge regression just split weight
# arbitrarily among them (near-perfect train AND test fit, even/odd-frame
# split leaks via adjacency, no real signal). v2 instead uses each cell's own
# RAW, un-spliced existence window as a separate feature (1332 candidates,
# naturally non-redundant since siblings occupy different positions once
# divided), and validates with a genuine block-level train/test split
# (train on epochs 1,3,5; test on held-out epochs 2,4,6) instead of
# alternating frames.

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

const MIN_DIST = 10f0
norm3(v) = sqrt(sum(abs2, v))

println("Building raw (un-spliced) direction tracks for $(length(grouped)) distinct cell names...")
dir_tracks = Dict{String,Dict{Int,Vector{Float64}}}()
for (name, rows) in grouped
    d = Dict{Int,Vector{Float64}}()
    for (t, p) in rows
        v = p - centroid_of[t]
        norm3(v) < MIN_DIST && continue
        d[t] = Float64.(normalize(v))
    end
    length(d) >= 10 && (dir_tracks[name] = d)  # need enough frames to be a useful feature
end
println("$(length(dir_tracks)) names retained (>=10 usable frames each).")

const NEPOCHS = 6
epoch_bounds = round.(Int, range(minimum(times), maximum(times), length=NEPOCHS + 1))
epochs = [epoch_bounds[i]:epoch_bounds[i+1] for i in 1:NEPOCHS]
train_times = reduce(vcat, [collect(epochs[i]) for i in 1:2:NEPOCHS])
test_times = reduce(vcat, [collect(epochs[i]) for i in 2:2:NEPOCHS])
train_times = unique(train_times); test_times = unique(test_times)
println("train frames: $(length(train_times)) (epochs 1,3,5)   test frames: $(length(test_times)) (epochs 2,4,6)")

function fit_axis(axis_field::Symbol, lambda::Float64; fit_times=train_times)
    names = collect(keys(dir_tracks))
    n = length(names)

    Xtr = zeros(3 * length(fit_times), n)
    ytr = zeros(3 * length(fit_times))
    for (ti, t) in enumerate(fit_times)
        target = getfield(axes_of[t], axis_field)
        for k in 1:3
            ytr[3*(ti-1)+k] = target[k]
        end
        for (ci, name) in enumerate(names)
            d = dir_tracks[name]
            haskey(d, t) || continue
            for k in 1:3
                Xtr[3*(ti-1)+k, ci] = d[t][k]
            end
        end
    end

    w = (Xtr'Xtr + lambda * LinearAlgebraI(n)) \ (Xtr'ytr)

    function fit_quality(ts)
        sims = Float64[]
        for t in ts
            target = getfield(axes_of[t], axis_field)
            pred = zeros(3)
            for (ci, name) in enumerate(names)
                d = dir_tracks[name]
                haskey(d, t) || continue
                pred .+= w[ci] .* d[t]
            end
            nrm = norm3(pred)
            nrm < 1e-6 && continue
            push!(sims, dot(pred ./ nrm, Vector(target)))
        end
        isempty(sims) ? NaN : mean(sims)
    end

    return names, w, fit_quality(train_times), fit_quality(test_times)
end

const LAMBDAS = [10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0]

# Final production weights per axis (refit on ALL 361 frames, at the lambda
# that generalized best in the train/test validation above) -- the
# validation split's job is done once we've picked lambda; the actual
# weights used downstream (CSV export, movie) should use every frame.
production_weights = Dict{Symbol,Any}()

for axis_field in (:dv, :lr)
    println()
    println("="^70)
    println("Axis target: Approach 1's $(axis_field)(t)  (single global fit, all 361 frames)")
    println("="^70)
    best = nothing
    for lambda in LAMBDAS
        names, w, train_fit, test_fit = fit_axis(axis_field, lambda)
        println("  lambda=$lambda: train_fit=$(round(train_fit,digits=3)) test_fit=$(round(test_fit,digits=3))")
        if best === nothing || test_fit > best.test_fit
            best = (; lambda, names, w, train_fit, test_fit)
        end
    end
    println("Best: lambda=$(best.lambda) train_fit=$(round(best.train_fit,digits=3)) test_fit=$(round(best.test_fit,digits=3))")
    order = sortperm(abs.(best.w), rev=true)
    println("Top 15 weighted cells:")
    for i in order[1:15]
        println("    $(best.names[i]): w=$(round(best.w[i],digits=4))")
    end

    names_full, w_full, train_fit_full, _ = fit_axis(axis_field, best.lambda; fit_times=times)
    println("Refit on all $(length(times)) frames at lambda=$(best.lambda): fit=$(round(train_fit_full,digits=3))")
    production_weights[axis_field] = (names=names_full, w=w_full)
end

# --- CSV export ---
# Both axes were fit over the same candidate set (`dir_tracks`, unchanged
# between calls), so `names` is in the same order for both -- safe to zip
# into one table.
dv_names, dv_w = production_weights[:dv].names, production_weights[:dv].w
lr_names, lr_w = production_weights[:lr].names, production_weights[:lr].w
@assert dv_names == lr_names "candidate name order differs between dv/lr fits"

csv_path = joinpath(@__DIR__, "..", "pretwitch_axis_regression_weights.csv")
open(csv_path, "w") do io
    println(io, "name,n_frames,dv_weight,lr_weight")
    order = sortperm(abs.(dv_w) .+ abs.(lr_w), rev=true)
    for i in order
        name = dv_names[i]
        n_frames = length(dir_tracks[name])
        println(io, "$name,$n_frames,$(dv_w[i]),$(lr_w[i])")
    end
end
println()
println("wrote $csv_path")
