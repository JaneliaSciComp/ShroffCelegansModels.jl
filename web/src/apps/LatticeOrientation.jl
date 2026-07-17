"""
    LatticeOrientation

Submodule for the lattice orientation QC page (deployment container
`lattice-orientation`, port 9401). Pure Bonito DOM — no Makie. Reads
`lattice_orientation.h5` (written by the check-lattice-orientation CronJob,
`web/scripts/check_lattice_orientation.jl`) and renders, per dataset, every
resolvable named check (`hyp7_Cpaaaa` plus the survey-derived candidate
cells) as a nested, highlighted pass/fail tree with sign, magnitude, and —
for inconsistent timepoints — a deep link into `/fix_annotation_ap_axis/`
for that exact dataset/cell/timepoint.

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

# A magnitude below this fraction of a dataset's own median magnitude (for
# that check) is flagged as a weak/low-confidence signal, independent of
# whether its sign happens to agree with the dataset's representative sign.
# Self-calibrated per (dataset, check) — no hardcoded absolute distance.
const LOW_CONFIDENCE_RATIO = 0.25

# Same URL scheme as `ShroffCelegansModels.get_fix_url` (src/demo_averaging/
# zscore_analysis.jl), reimplemented locally rather than reused: that
# function is defined via a runtime `Base.include` into the ZscoreAnalysis
# module (world-age hazard — see notes/), whereas this app must stay fully
# precompilable.
fix_ap_axis_url(group::AbstractString, group_idx::Integer, annotation::AbstractString, timepoint::Integer) =
    "https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/fix_annotation_ap_axis/$group/$group_idx?annotation=$annotation&timepoint=$timepoint"

struct AnnotationCheck
    annotation_name::String
    axis::Symbol
    signs::Vector{Float64}
    magnitudes::Vector{Float64}
    magnitude_median::Float64
    representative_sign::Float64
    reference_sign::Float64
    matches_reference::Bool
    mismatched_timepoint_count::Int
    link_eligible::Bool
end

struct DatasetOrientation
    index::Int
    path::String
    cell_key_name::String
    start::Int
    stop::Int
    outliers::Vector{Int}
    checks::Vector{AnnotationCheck}
end

function read_lattice_orientation(path::AbstractString)
    h5open(path, "r") do f
        root = f["lattice_orientation"]
        groups = Dict{String, Vector{DatasetOrientation}}()
        for group_name in keys(root)
            g = root[group_name]
            indices = sort(parse.(Int, collect(keys(g))))
            entries = map(indices) do i
                dsg = g[string(i)]
                A = attrs(dsg)
                check_count = Int(A["check_count"])
                checks = map(1:check_count) do c
                    cg = dsg[string("check_", c)]
                    CA = attrs(cg)
                    AnnotationCheck(
                        string(CA["annotation_name"]),
                        Symbol(CA["axis"]),
                        Float64.(read(cg["sign"])),
                        Float64.(CA["magnitude"]),
                        Float64(CA["magnitude_median"]),
                        Float64(CA["representative_sign"]),
                        Float64(CA["reference_sign"]),
                        Bool(CA["matches_reference"]),
                        Int(CA["mismatched_timepoint_count"]),
                        Bool(CA["link_eligible"]),
                    )
                end
                DatasetOrientation(
                    i,
                    string(A["path"]),
                    string(A["cell_key.name"]),
                    Int(A["cell_key.start"]),
                    Int(A["cell_key.end"]),
                    Int.(A["cell_key.outliers"]),
                    checks,
                )
            end
            groups[group_name] = entries
        end
        return groups
    end
end

format_sign(axis::Symbol, s::Real) = isnan(s) ? "no annotation" :
    axis === :dv ? (s > 0 ? "+1 (dorsal)" : "-1 (ventral)") : (s > 0 ? "+1 (right)" : "-1 (left)")
format_magnitude(m::Real) = isnan(m) ? "n/a" : string(round(m; digits=2))
axis_label(axis::Symbol) = axis === :dv ? "DV" : "LR"

const HIGHLIGHT_CLASS = "shroff-highlight"

# One row per distinct (annotation_name, axis) check, aggregated across every
# dataset that resolves it — lets the page header summarize all ~30 checks
# without repeating each dataset's full detail.
function check_overview(groups::Dict{String, Vector{DatasetOrientation}})
    seen = Dict{Tuple{String,Symbol}, @NamedTuple{reference_sign::Float64, resolved::Int, mismatched::Int}}()
    for datasets in values(groups), d in datasets, c in d.checks
        key = (c.annotation_name, c.axis)
        prev = get(seen, key, (reference_sign = c.reference_sign, resolved = 0, mismatched = 0))
        seen[key] = (
            reference_sign = c.reference_sign,
            resolved = prev.resolved + 1,
            mismatched = prev.mismatched + (c.matches_reference ? 0 : 1),
        )
    end
    return seen
end

function render_summary(groups::Dict{String, Vector{DatasetOrientation}}, source_mtime::Float64)
    overview = check_overview(groups)
    DOM.main(
        theme_assets()...,
        DOM.h2("Lattice orientation"),
        DOM.p(
            "Source: ", DOM.code(lattice_orientation_path()),
            " (file mtime: ", isnan(source_mtime) ? "" : string(source_mtime), ")",
        ),
        DOM.h3("Checks"),
        DOM.ul(map(sort(collect(keys(overview)); by = k -> (string(k[2]), k[1]))) do key
            name, axis = key
            o = overview[key]
            DOM.li(
                DOM.code(name), " (", axis_label(axis), ") — reference sign: ",
                DOM.strong(format_sign(axis, o.reference_sign)),
                " — ", string(o.resolved - o.mismatched), "/", string(o.resolved),
                " datasets agree",
            )
        end),
        map(sort(collect(keys(groups)))) do group
            DOM.section(
                DOM.h3(group),
                render_group(group, groups[group]),
            )
        end...,
        ; class="shroff-mtime container",
    )
end

function render_group(group::AbstractString, datasets::Vector{DatasetOrientation})
    DOM.ul(map(datasets) do d
        any_fail = any(c -> !c.matches_reference, d.checks)
        total_mismatched_timepoints = sum(c -> c.mismatched_timepoint_count, d.checks; init=0)
        DOM.li(DOM.details(
            DOM.summary(
                "[", string(d.index), "] ",
                DOM.code(d.cell_key_name),
                " — ", string(length(d.checks)), " check(s) resolved, ",
                string(total_mismatched_timepoints), " timepoint mismatch(es)",
                any_fail ? " (check for swap)" : "",
                ; class = any_fail ? HIGHLIGHT_CLASS : "",
            ),
            DOM.div("path: ", DOM.code(d.path)),
            DOM.div(
                "timepoints: ", string(d.start), "–", string(d.stop),
                ", outliers: ", isempty(d.outliers) ? "none" : join(d.outliers, ", "),
            ),
            isempty(d.checks) ?
                DOM.p("No checked cells were resolvable for this dataset.") :
                DOM.ul(map(c -> render_check(group, d, c), d.checks)),
        ))
    end...)
end

function render_check(group::AbstractString, d::DatasetOrientation, c::AnnotationCheck)
    n_timepoints = length(c.signs)
    status_label = c.matches_reference ? " (PASS)" : " (FAIL — check for swap)"
    DOM.li(DOM.details(
        DOM.summary(
            DOM.code(c.annotation_name), " (", axis_label(c.axis), ")",
            " — representative sign: ", format_sign(c.axis, c.representative_sign),
            status_label,
            " — ", string(c.mismatched_timepoint_count), "/", string(n_timepoints), " timepoints mismatched",
            " — magnitude (median): ", format_magnitude(c.magnitude_median),
            ; class = c.matches_reference ? "" : HIGHLIGHT_CLASS,
        ),
        DOM.ul(map(eachindex(c.signs)) do i
            tp = d.start + i - 1
            s = c.signs[i]
            mag = c.magnitudes[i]
            inconsistent = !isnan(s) && s != c.representative_sign
            weak = !isnan(mag) && !isnan(c.magnitude_median) && c.magnitude_median > 0 &&
                mag < LOW_CONFIDENCE_RATIO * c.magnitude_median
            tp_class = (inconsistent || weak) ? HIGHLIGHT_CLASS : ""
            if isnan(s)
                DOM.li("t=", string(tp), ": ", tp in d.outliers ? "outlier" : "no annotation"; class=tp_class)
            else
                label = join(
                    filter(!isempty, [
                        format_sign(c.axis, s) * " (mag " * format_magnitude(mag) * ")",
                        inconsistent ? "inconsistent with dataset" : "",
                        weak ? "weak signal" : "",
                    ]),
                    " — ",
                )
                fix_link = (inconsistent && c.link_eligible) ?
                    DOM.a(" [fix]"; href=fix_ap_axis_url(group, d.index, c.annotation_name, tp), target="_blank") :
                    ""
                DOM.li("t=", string(tp), ": ", label, fix_link; class=tp_class)
            end
        end),
    ))
end

function render_missing(path::AbstractString)
    DOM.div(
        DOM.h2("Lattice orientation"),
        DOM.p("No data yet. Expected ", DOM.code(path), " but it does not exist."),
        DOM.p("This file is generated by the check-lattice-orientation CronJob."),
    )
end

function web_lattice_orientation()
    server = Server(
        "0.0.0.0", PORT;
        proxy_url="https://$(get(ENV, "SHROFF_HOST", "shroff-data.int.janelia.org"))/lattice_orientation/",
    )
    route!(server, "/" => App(; title="Shroff C. elegans lattice orientation") do
        path = lattice_orientation_path()
        if !isfile(path)
            return render_missing(path)
        end
        groups = read_lattice_orientation(path)
        return render_summary(groups, Float64(mtime(path)))
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
    mk_check(name, axis, rep, matches; magnitudes=Float64[3.5, NaN, 3.1], link_eligible=true) = AnnotationCheck(
        name, axis,
        Float64[rep, NaN, matches ? rep : -rep],
        magnitudes, 3.3, rep, rep, matches, matches ? 0 : 1, link_eligible,
    )
    mk(i, rep, matches) = DatasetOrientation(
        i, "/nearline/shroff/example/Pos$i/RegB", "cellkey$i", 1, 3, Int[],
        [mk_check("hyp7_Cpaaaa", :dv, rep, matches), mk_check("AVDL", :lr, -1.0, true; link_eligible=true)],
    )
    # A dataset with no resolvable checks at all.
    no_data = DatasetOrientation(
        2, "/nearline/shroff/example/Pos2/RegB", "cellkey2", 1, 3, Int[], AnnotationCheck[],
    )
    Dict{String, Vector{DatasetOrientation}}(
        "RW10000" => [mk(0, 1.0, true), mk(1, -1.0, false), no_data],
    )
end

if @load_preference("precompile_workload", true)
@setup_workload begin
    groups = _synthetic_groups()
    @compile_workload begin
        try
            app = App(() -> render_summary(groups, 1.7e9))
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
