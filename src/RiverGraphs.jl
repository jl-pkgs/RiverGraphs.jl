module RiverGraphs


using Graphs, Parameters, DataFrames
import SpatRasters: st_dims
import RTableTools: cbind

include("IO.jl")
include("RiverGraph.jl")
include("fillnodata_upstream.jl")

include("flow_path.jl")
include("subbasins.jl")
include("stream.jl")
include("subdomains.jl")

include("utils.jl")
include("sf.jl")


export RiverGraph
export active_indices, reverse_index, pcr_dir,
  graph_flow,
  topological_sort_by_dfs,
  stream_order, stream_link, stream_network,
  fillnodata_upstream, fillnodata_upriver

const path_flowdir_GuanShan = abspath("$(@__DIR__)/../data/GuanShan_flwdir.tif")


export path_flowdir_GuanShan
export nv # from Graphs


end # module RiverGraphs
