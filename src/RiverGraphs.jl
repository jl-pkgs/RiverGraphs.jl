module RiverGraphs


using Graphs, Parameters, DataFrames
using Printf
import SpatialRasterLite: st_dims
import RTableTools: cbind, fwrite, fread

include("IO.jl")
include("RiverGraph.jl")
include("constrained_flowdir.jl")
include("priority_flood.jl")
include("fillnodata_upbasin.jl")
include("fillnodata_upriver.jl")

include("flow_path.jl")
include("st_watershed.jl")
include("stream_order.jl")
include("stream_link.jl")
include("st_stream_network.jl")
include("subdomains.jl")

include("utils.jl")
include("sf.jl")


"""
    read_river_lines(source; kwargs...)

Read river centreline geometries as vectors of `(x, y)` coordinate tuples.
Methods for GIS vector files are provided by the optional ArchGDAL extension;
load ArchGDAL before calling this function on a file path.
"""
function read_river_lines end


export RiverGraph
export active_indices, reverse_index, pcr_dir,
  graph_flow,
  topological_sort_kahn,
  stream_order, stream_link, stream_network,
  fillnodata_upbasin, fillnodata_upriver
export read_river_lines

const path_flowdir_GuanShan = abspath("$(@__DIR__)/../data/GuanShan_flwdir.tif")


export path_flowdir_GuanShan
export SimpleDiGraph, nv # from Graphs
export DataFrame, nrow


end # module RiverGraphs
