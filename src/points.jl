module Points

using ..Types: AbstractCelegansModel, lattice_knots, transverse_splines, names
export points, lattice, cross_sections, CrossSections

function points(model::AbstractCelegansModel)
    splines = transverse_splines(model)
    r = LinRange(0, 1, length(model))
    pts = map(splines) do spline
        spline.(r)
    end
    hcat(pts...)
end

function lattice(model::AbstractCelegansModel)
    splines = transverse_splines(model)
    r = lattice_knots(model)
    pts = map(splines) do spline
        spline.(r)
    end
    hcat(pts...)
end

struct CrossSections{P <: Vector{<: Pair}}
    sections::P
end

function cross_sections(model::AbstractCelegansModel)
    _names = names(model)[1:2:end]
    CrossSections(_names .=> copy.(eachrow(lattice(model))))
end


end