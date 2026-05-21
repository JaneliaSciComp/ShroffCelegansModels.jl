using Bonito
using HDF5: h5open, attrs
using Dates: unix2datetime, format as date_format

const PORT = 9400

modified_times_path() = joinpath(
    get(ENV, "MODIFIED_TIMES_DIR", "/data/annotations/modified_times"),
    "modified_times.h5",
)

struct DatasetEntry
    index::Int
    path::String
    cell_key_name::String
    start::Int
    stop::Int
    outliers::Vector{Int}
    mtimes::Vector{Float64}
end

function read_modified_times(path::AbstractString)
    h5open(path, "r") do f
        groups = Dict{String, Vector{DatasetEntry}}()
        for group_name in keys(f)
            g = f[group_name]
            indices = sort(parse.(Int, collect(keys(g))))
            entries = map(indices) do i
                ds = g[string(i)]
                A = attrs(ds)
                DatasetEntry(
                    i,
                    string(A["path"]),
                    string(A["cell_key.name"]),
                    Int(A["cell_key.start"]),
                    Int(A["cell_key.end"]),
                    Int.(A["cell_key.outliers"]),
                    Float64.(read(ds)),
                )
            end
            groups[group_name] = entries
        end
        return groups
    end
end

format_unix(u::Real) = isnan(u) ? "" : date_format(unix2datetime(u), "yyyy-mm-dd HH:MM:SS")

dataset_max(d::DatasetEntry) = (v = filter(!isnan, d.mtimes); isempty(v) ? -Inf : maximum(v))

function last_modified(mtimes::Vector{Float64})
    valid = filter(!isnan, mtimes)
    isempty(valid) ? "n/a" : format_unix(maximum(valid))
end

# Indices whose value equals the (finite) maximum of `xs`. Empty if all NaN.
function argmax_finite(xs::Vector{Float64})
    valid = filter(!isnan, xs)
    isempty(valid) && return Int[]
    m = maximum(valid)
    return findall(x -> !isnan(x) && x == m, xs)
end

const HIGHLIGHT_STYLE = "background:#fffbcc;font-weight:bold"

function render_groups(groups::Dict{String, Vector{DatasetEntry}}, source_mtime::Float64)
    DOM.div(
        DOM.h2("Integrated annotation modified times"),
        DOM.p(
            "Source: ", DOM.code(modified_times_path()),
            " (file mtime: ", format_unix(source_mtime), ")",
        ),
        let
            sorted_names = sort(collect(keys(groups)))
            group_maxes = Dict(name => let
                ms = dataset_max.(groups[name])
                isempty(ms) ? -Inf : maximum(ms)
            end for name in sorted_names)
            global_max = isempty(group_maxes) ? -Inf : maximum(values(group_maxes))
            map(sorted_names) do group_name
                datasets = groups[group_name]
                ds_maxes = dataset_max.(datasets)
                group_max = group_maxes[group_name]
                hot_datasets = Set(findall(==(group_max), ds_maxes))
                group_last = group_max == -Inf ? "n/a" : format_unix(group_max)
                group_summary_style = group_max == global_max ? HIGHLIGHT_STYLE : ""
                DOM.details(
                    DOM.summary(
                        DOM.strong(group_name),
                        " — ", string(length(datasets)), " datasets",
                        " — last modified: ", group_last,
                        ; style=group_summary_style,
                    ),
                    DOM.ul(map(enumerate(datasets)) do (ds_i, d)
                        hot_tps = Set(argmax_finite(d.mtimes))
                        dataset_summary_style = ds_i in hot_datasets ? HIGHLIGHT_STYLE : ""
                        DOM.li(DOM.details(
                            DOM.summary(
                                "[", string(d.index), "] ",
                                DOM.code(d.cell_key_name),
                                " — last modified: ", last_modified(d.mtimes),
                                ; style=dataset_summary_style,
                            ),
                            DOM.div("path: ", DOM.code(d.path)),
                            DOM.div(
                                "timepoints: ", string(d.start), "–", string(d.stop),
                                ", outliers: ", isempty(d.outliers) ? "none" : join(d.outliers, ", "),
                            ),
                            DOM.ul(map(eachindex(d.mtimes)) do i
                                tp = d.start + i - 1
                                u = d.mtimes[i]
                                label = if isnan(u)
                                    tp in d.outliers ? "outlier" : "missing"
                                else
                                    format_unix(u)
                                end
                                tp_style = i in hot_tps ? HIGHLIGHT_STYLE : ""
                                DOM.li("t=", string(tp), ": ", label; style=tp_style)
                            end),
                        ))
                    end),
                )
            end
        end,
    )
end

function render_missing(path::AbstractString)
    DOM.div(
        DOM.h2("Integrated annotation modified times"),
        DOM.p("No data yet. Expected ", DOM.code(path), " but it does not exist."),
        DOM.p("This file is generated by the save-modified-times CronJob."),
    )
end

function web_modified_times()
    server = Server(
        "0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/modified_times/",
    )
    route!(server, "/" => App(; title="Shroff C. elegans modified times") do
        path = modified_times_path()
        if !isfile(path)
            return render_missing(path)
        end
        groups = read_modified_times(path)
        return render_groups(groups, Float64(mtime(path)))
    end)
    return server
end

function web_main()
    server = web_modified_times()
    @info "Listening" port=PORT
    return server
end

if abspath(PROGRAM_FILE) == @__FILE__
    wait(web_main())
end
