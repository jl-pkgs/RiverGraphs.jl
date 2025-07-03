module RiverGraphs


using Graphs, Parameters, DataFrames
import SpatRasters: st_dims
import RTableTools: cbind

include("IO.jl")
include("RiverGraph.jl")
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


export RiverGraph
export active_indices, reverse_index, pcr_dir,
  graph_flow,
  topological_sort_by_dfs,
  stream_order, stream_link, stream_network,
  fillnodata_upbasin, fillnodata_upriver

const path_flowdir_GuanShan = abspath("$(@__DIR__)/../data/GuanShan_flwdir.tif")


export path_flowdir_GuanShan
export nv # from Graphs


end # module RiverGraphs
