"""
    LatticeOrientation

Submodule for the lattice LR-orientation QC page (deployment container
`lattice-orientation`, port 9401). Pure Bonito DOM — no Makie. Reads
`lattice_orientation.h5` (written by the check-lattice-orientation CronJob)
and renders a nested, highlighted view of per-dataset/per-timepoint
Cpaaaa-vs-seam-plane orientation signs, magnitudes, and LR cross-checks.

The precompile workload renders a tiny synthetic tree through
`Bonito.export_static`, caching the DOM-serialize path with no disk dependency
(the only asset is an in-tree `web/static/style.css` plus a CDN URL).
"""
module LatticeOrientation

using Bonito
using Bonito: Asset
using HDF5: h5open, attrs
using PrecompileTools: @setup_workload, @compile_workload
using Preferences: @load_preference

const PORT = 9401

# --- theme assets (moved from web/scripts/web_theme.jl) ------------------------
# pico.css v2 (classless, dark/light auto) from CDN + a small local override.
const PICO_CSS = Asset("https://cdn.jsdelivr.net/npm/@picocss/pico@2/css/pico.min.css")
const SHROFF_CSS = Asset(joinpath(@__DIR__, "..", "..", "static", "style.css"))
theme_assets() = (PICO_CSS, SHROFF_CSS)

lattice_orientation_path() = joinpath(
    get(ENV, "LATTICE_ORIENTATION_DIR", "/data/annotations/lattice_orientation"),
    "lattice_orientation.h5",
)

# A DV (or LR) magnitude below this fraction of a dataset's own median
# magnitude is flagged as a weak/low-confidence signal, independent of
# whether its sign happens to agree with the dataset's representative sign.
# Self-calibrated per dataset — no hardcoded absolute distance.
const LOW_CONFIDENCE_RATIO = 0.25

struct DatasetOrientation
    index::Int
    path::String
    cell_key_name::String
    start::Int
    stop::Int
    outliers::Vector{Int}
    dv_signs::Vector{Float64}
    dv_magnitudes::Vector{Float64}
    dv_magnitude_median::Float64
    lr_signs::Vector{Float64}
    lr_magnitudes::Vector{Float64}
    lr_representative_sign::Float64
    lr_magnitude_median::Float64
    cpaaaa_key::String
    representative_sign::Float64
    matches_reference::Bool
    mismatched_timepoint_count::Int
end

function read_lattice_orientation(path::AbstractString)
    h5open(path, "r") do f
        reference_sign = Float64(attrs(f)["reference_sign"])
        root = f["lattice_orientation"]
        groups = Dict{String, Vector{DatasetOrientation}}()
        for group_name in keys(root)
            g = root[group_name]
            indices = sort(parse.(Int, collect(keys(g))))
            entries = map(indices) do i
                ds = g[string(i)]
                A = attrs(ds)
                DatasetOrientation(
                    i,
                    string(A["path"]),
                    string(A["cell_key.name"]),
                    Int(A["cell_key.start"]),
                    Int(A["cell_key.end"]),
                    Int.(A["cell_key.outliers"]),
                    Float64.(read(ds)),
                    Float64.(A["dv_magnitude"]),
                    Float64(A["dv_magnitude_median"]),
                    Float64.(A["lr_sign"]),
                    Float64.(A["lr_magnitude"]),
                    Float64(A["lr_representative_sign"]),
                    Float64(A["lr_magnitude_median"]),
                    string(A["cpaaaa_key"]),
                    Float64(A["representative_sign"]),
                    Bool(A["matches_reference"]),
                    Int(A["mismatched_timepoint_count"]),
                )
            end
            groups[group_name] = entries
        end
        return reference_sign, groups
    end
end

format_sign(s::Real) = isnan(s) ? "no annotation" : (s > 0 ? "+1 (dorsal)" : "-1 (ventral)")
format_lr_sign(s::Real) = isnan(s) ? "n/a" : (s > 0 ? "+1 (right)" : "-1 (left)")
format_magnitude(m::Real) = isnan(m) ? "n/a" : string(round(m; digits=2))

const HIGHLIGHT_CLASS = "shroff-highlight"

function render_summary(reference_sign::Float64, groups::Dict{String, Vector{DatasetOrientation}}, source_mtime::Float64)
    all_entries = collect(Iterators.flatten(values(groups)))
    mismatched = count(d -> !d.matches_reference, all_entries)
    total_mismatched_timepoints = sum(d -> d.mismatched_timepoint_count, all_entries; init=0)
    DOM.main(
        theme_assets()...,
        DOM.h2("Lattice LR orientation"),
        DOM.p(
            "Source: ", DOM.code(lattice_orientation_path()),
            " (file mtime: ", isnan(source_mtime) ? "" : string(source_mtime), ")",
        ),
        DOM.p(
            "Reference (dorsal) sign: ", DOM.strong(format_sign(reference_sign)),
            " — ", string(length(all_entries) - mismatched), "/", string(length(all_entries)),
            " datasets agree; ", string(total_mismatched_timepoints),
            " individual timepoint(s) disagree with their own dataset across all datasets",
        ),
        map(sort(collect(keys(groups)))) do group
            DOM.section(
                DOM.h3(group),
                render_group(groups[group], reference_sign),
            )
        end...,
        ; class="shroff-mtime container",
    )
end

function render_group(datasets::Vector{DatasetOrientation}, reference_sign::Float64)
    DOM.ul(map(datasets) do d
        dataset_class = d.matches_reference ? "" : HIGHLIGHT_CLASS
        n_timepoints = length(d.dv_signs)
        DOM.li(DOM.details(
            DOM.summary(
                "[", string(d.index), "] ",
                DOM.code(d.cell_key_name),
                " — representative sign: ", format_sign(d.representative_sign),
                d.matches_reference ? " (PASS)" : " (FAIL — check for LR swap)",
                " — ", string(d.mismatched_timepoint_count), "/", string(n_timepoints),
                " timepoints mismatched",
                ; class=dataset_class,
            ),
            DOM.div("path: ", DOM.code(d.path)),
            DOM.div("Cpaaaa annotation key: ", DOM.code(isempty(d.cpaaaa_key) ? "not found" : d.cpaaaa_key)),
            DOM.div(
                "timepoints: ", string(d.start), "–", string(d.stop),
                ", outliers: ", isempty(d.outliers) ? "none" : join(d.outliers, ", "),
            ),
            DOM.div(
                "DV magnitude (median): ", format_magnitude(d.dv_magnitude_median),
                " — LR: representative sign ", format_lr_sign(d.lr_representative_sign),
                ", magnitude (median) ", format_magnitude(d.lr_magnitude_median),
            ),
            DOM.ul(map(eachindex(d.dv_signs)) do i
                tp = d.start + i - 1
                s = d.dv_signs[i]
                mag = d.dv_magnitudes[i]
                lr_s = d.lr_signs[i]
                lr_mag = d.lr_magnitudes[i]
                inconsistent = !isnan(s) && !isnan(d.representative_sign) && s != d.representative_sign
                weak = !isnan(mag) && !isnan(d.dv_magnitude_median) && d.dv_magnitude_median > 0 &&
                    mag < LOW_CONFIDENCE_RATIO * d.dv_magnitude_median
                tp_class = (inconsistent || weak) ? HIGHLIGHT_CLASS : ""
                label = if isnan(s)
                    tp in d.outliers ? "outlier" : "no Cpaaaa annotation"
                else
                    join(
                        filter(!isempty, [
                            format_sign(s) * " (mag " * format_magnitude(mag) * ")",
                            inconsistent ? "inconsistent with dataset" : "",
                            weak ? "weak signal" : "",
                            "LR " * format_lr_sign(lr_s) * " (mag " * format_magnitude(lr_mag) * ")",
                        ]),
                        " — ",
                    )
                end
                DOM.li("t=", string(tp), ": ", label; class=tp_class)
            end),
        ))
    end...)
end

function render_missing(path::AbstractString)
    DOM.div(
        DOM.h2("Lattice LR orientation"),
        DOM.p("No data yet. Expected ", DOM.code(path), " but it does not exist."),
        DOM.p("This file is generated by the check-lattice-orientation CronJob."),
    )
end

function web_lattice_orientation()
    server = Server(
        "0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/lattice_orientation/",
    )
    route!(server, "/" => App(; title="Shroff C. elegans lattice LR orientation") do
        path = lattice_orientation_path()
        if !isfile(path)
            return render_missing(path)
        end
        reference_sign, groups = read_lattice_orientation(path)
        return render_summary(reference_sign, groups, Float64(mtime(path)))
    end)
    return server
end

function main()
    server = web_lattice_orientation()
    @info "Listening" port=PORT
    if isinteractive()
        println("Press enter to quit")
        readline()
    else
        wait(server)
    end
end

# --- precompile workload -------------------------------------------------------

function _synthetic_groups()
    mk(i, rep, matches) = DatasetOrientation(
        i, "/nearline/shroff/example/Pos$i/RegB", "cellkey$i", 1, 3,
        Int[],
        Float64[1.0, NaN, matches ? 1.0 : -1.0],
        Float64[3.5, NaN, 3.1],
        3.3,
        Float64[-1.0, NaN, -1.0],
        Float64[0.4, NaN, 0.5],
        -1.0,
        0.45,
        "hyp7_Cpaaaa", rep, matches, matches ? 0 : 1,
    )
    Dict{String, Vector{DatasetOrientation}}(
        "RW10000" => [mk(0, 1.0, true), mk(1, -1.0, false)],
    )
end

if @load_preference("precompile_workload", true)
@setup_workload begin
    groups = _synthetic_groups()
    @compile_workload begin
        try
            app = App(() -> render_summary(1.0, groups, 1.7e9))
            mktempdir() do dir
                export_static(joinpath(dir, "lattice_orientation.html"), app)
            end
        catch err
            @debug "LatticeOrientation precompile workload skipped" exception = (err, catch_backtrace())
        end
    end
end
end  # precompile_workload preference guard

end # module LatticeOrientation
