module Datasets
    using JSON3: JSON3

    export Dataset, CellKey, NormalizedDataset

    abstract type AbstractDataset end

    const CELL_KEY_FILE_NAME = "cell_key.json"

    struct Dataset{JO <: JSON3.Object} <: AbstractDataset
        path::String
        cell_key::JO
    end
    function Dataset(path::String)
        path = abspath(path)
        isdir(path) || SystemError("Dataset path, \"$path\", does not exist as a directory.")
        cell_key_path = joinpath(dirname(path), CELL_KEY_FILE_NAME)
        isfile(cell_key_path) || SystemError("Cell Key path, $cell_key_path, does not exist as a file.")
        return _Dataset(path, cell_key_path)
    end
    function _Dataset(path::String, cell_key_path::String)
        return Dataset(path, JSON3.read(cell_key_path))
    end

    struct CellKey
        name::String
        start::Int
        stop::Int # end
        mapping::Dict{Symbol,String}
        outliers::Vector{Int}
    end
    const _cell_key_keys = (:end, :mapping, :name, :outliers, :start)
    function CellKey(jo::JSON3.Object)
        _cell_key_keys ⊆ keys(jo) || return ArgumentError("JSON3.Object must contain the following keys: $_cell_key_keys")
        CellKey(jo.name, jo.start, jo.end, jo.mapping, jo.outliers)
    end
    function Base.getproperty(ck::CellKey, s::Symbol)
        if s == :end
            s = :stop
        end
        return Base.getfield(ck, s)
    end
    function Base.range(ck::CellKey)
        return range(ck.start, ck.stop)
    end
    
    struct NormalizedDataset <: AbstractDataset
        path::String
        cell_key::CellKey
    end
    function NormalizedDataset(path::String)
        ds = Dataset(path)
        return NormalizedDataset(ds)
    end
    function NormalizedDataset(ds::Dataset)
        return NormalizedDataset(ds.path, CellKey(ds.cell_key))
    end
end