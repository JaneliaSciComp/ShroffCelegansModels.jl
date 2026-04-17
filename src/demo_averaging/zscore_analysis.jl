using HDF5
using StatsBase: zscore
using DataFrames
using CSV

function zscore_analysis()
    df = DataFrame(dataset=String[], embryo=String[], annotation=String[], zscore=Float64[])
    h5open("embryos_371_2026_04_10.h5") do h5f
        for k in keys(h5f)
            for e in keys(h5f[k])
                M = h5f[k][e][3, :, :]
                zscores = zscore(sqrt.(sum(diff(M; dims=1) .^ 2; dims=1)))
                pairs = Dict(attrs(h5f[k][e])["annotation_names"] .=> zscores')
                deviation = filter(k -> pairs[k] > 2, keys(pairs))
                for d in deviation
                    println(k, ", ", e, ", ", d, ", ", pairs[d])
                    push!(df, (k, e, d, pairs[d]))
                end
            end
        end
    end
    return df
end