function Base.show(io::IO, m::MIME"text/plain", model::T) where T <: AbstractCelegansModel
    println(io, T)
    for f in fieldnames(T)
        println(io, "    ", f)
    end
end