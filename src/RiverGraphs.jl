module RiverGraphs


using Graphs, Parameters, DataFrames
import Ipaper: read_flowdir
import Ipaper.sf: st_dims
import RTableTools: cbind

include("RiverGraph.jl")
include("flow_path.jl")
include("subbasins.jl")
include("stream.jl")
include("utils.jl")


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
