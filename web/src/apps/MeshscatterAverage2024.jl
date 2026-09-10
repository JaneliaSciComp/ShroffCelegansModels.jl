"""
    MeshscatterAverage2024

Web app submodule for the `meshscatter_average_edited_2024_10_24` viewer
(deployment container `meshscatter-2024-10-24`, port 8590).

Identical render to [`MeshscatterAverage`](@ref); it differs only in the source
HDF5 filename, port, and proxy path. It carries **no precompile workload** — the
render specializations are identical to MeshscatterAverage's and already cached.
"""
module MeshscatterAverage2024

using ..MeshscatterAverage: MeshscatterAverage

const DEFAULT_FILENAME =
    "edited_smoothed_average_annotations_r020_theta020_z030_with_seam_cells.h5"
const PORT = 8590
const PROXY_PATH = "meshscatter_average_edited_2024_10_24"

main() = MeshscatterAverage.main(; default_filename = DEFAULT_FILENAME, port = PORT, proxy_path = PROXY_PATH)

end # module MeshscatterAverage2024
