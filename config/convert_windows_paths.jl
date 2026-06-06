#!/usr/bin/env julia
# Converts Windows-style paths in a config JSON to Unix-style paths.
# Usage: julia convert_windows_paths.jl [input.json [output.json]]
#
# Default input:  config/windows/config_2026_03_19_v2.json
# Default output: config/linux/config_2026_03_19_v2.json
#
# Path mapping applied:
#   X:\shrofflab\  →  /nearline/shroff/shrofflab/
#   remaining backslashes  →  forward slashes

using JSON3

const WINDOWS_PREFIX = Regex("^[A-Za-z]:\\\\shrofflab\\\\", "i")
const UNIX_PREFIX    = "/nearline/shroff/shrofflab/"

function convert_path(s::AbstractString)
    s = replace(s, WINDOWS_PREFIX => UNIX_PREFIX)
    replace(s, '\\' => '/')
end

function convert_value(v)
    if v isa AbstractString
        return convert_path(v)
    elseif v isa AbstractDict
        return Dict(k => convert_value(val) for (k, val) in v)
    elseif v isa AbstractArray
        return [convert_value(el) for el in v]
    else
        return v
    end
end

function main()
    script_dir = dirname(abspath(@__FILE__))

    input_path = length(ARGS) >= 1 ? ARGS[1] :
        joinpath(script_dir, "windows", "config_2026_03_19_v2.json")

    default_output_dir = joinpath(script_dir, "linux")
    output_path = length(ARGS) >= 2 ? ARGS[2] :
        joinpath(default_output_dir, basename(input_path))

    isfile(input_path) || error("Input file not found: $input_path")

    config = open(input_path) do io
        JSON3.read(io, Dict{String,Any})
    end
    converted = convert_value(config)

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        JSON3.write(io, converted)
        println(io)           # trailing newline
    end

    println("Wrote converted config to: $output_path")
end

main()
