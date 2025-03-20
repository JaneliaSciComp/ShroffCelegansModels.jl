module Datasets
    using JSON3: JSON3
    using CSV, DataFrames

    export Dataset, CellKey, NormalizedDataset

    abstract type AbstractDataset end

    const CELL_KEY_FILE_NAME = "cell_key.json"
    const CELL_KEY_CSV_NAME = "CellKey.csv"

    struct Dataset{JO <: JSON3.Object} <: AbstractDataset
        path::String
        cell_key::JO
    end
    function Dataset(path::String)
        path = abspath(path)
        isdir(path) || error("Dataset path, \"$path\", does not exist as a directory.")
        cell_key_path = joinpath(dirname(path), CELL_KEY_FILE_NAME)
        # isfile(cell_key_path) || error("Cell Key path, $cell_key_path, does not exist as a file.")
        if isfile(cell_key_path)
            return _Dataset(path, cell_key_path)
        end
        cell_key_path = joinpath(dirname(path), CELL_KEY_CSV_NAME)
        return _Dataset(path, cell_key_path)
    end
    function _Dataset(path::String, cell_key_path::String)
        if endswith(cell_key_path, ".json")
            return Dataset(path, JSON3.read(cell_key_path))
        end
        if endswith(cell_key_path, ".csv")
            return NormalizedDataset(path, parse_cell_key_csv(cell_key_path))
        end
        error("Cell Key path, $cell_key_path, does not exist as a file.")
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
    function CellKey(df::DataFrame)
        name = df[1,1]
        start = df[2,1] isa AbstractString ? parse(Int, df[2,1]) : df[2,1]
        stop = df[2,2] isa AbstractString ? parse(Int, df[2,2]) : df[2,2]
        mapping = Dict{Symbol,String}(Symbol(k) => v for (k,v) in zip(df[4:end,1], df[4:end,2]))
        outliers = Int[]
        return CellKey(name, start, stop, mapping, outliers)
    end
    function Base.getproperty(ck::CellKey, s::Symbol)
        if s == :end
            s = :stop
        end
        return Base.getfield(ck, s)
    end
    function Base.range(ck::CellKey)
        r = range(ck.start, ck.stop)
        rd = setdiff(r, ck.outliers)
        real_start = minimum(rd)
        real_stop = maximum(rd)
        return range(real_start, real_stop)
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
    function NormalizedDataset(ds::NormalizedDataset)
        return NormalizedDataset(ds.path, ds.cell_key)
    end
    function parse_cell_key_csv(path::AbstractString)
        df = CSV.read(path, DataFrame, header=0)
        return CellKey(df)
    end
end
