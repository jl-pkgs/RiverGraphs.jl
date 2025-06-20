
using Graphs, Parameters
import Ipaper.sf: st_dims
import Ipaper: read_flowdir


"""
    graph_children(net, v)

```julia
julia> graph_children(net, 12)
Any[1, 2, [3, 4, 5, 7] => 8, [6, 9, 10] => 11] => 12
```
"""
function graph_children(net::SimpleDiGraph, v::Int)
  nodes_up = inneighbors(net, v) # 上游的站点
  isempty(nodes_up) && return v
  nodes = map(v -> graph_children(net, v), nodes_up) # 
  nodes => v
end

function Base.show(io::IO, net::SimpleDiGraph)
  v = find_outlet(net)
  println(graph_children(net, v))
end

find_outlet(net::SimpleDiGraph) = topological_sort_by_dfs(net)[end]

# 其他函数调用时，需要知道数据类型，因此此处保留了{FT}
@with_kw mutable struct RiverGraph{FT}
  data::AbstractArray{FT} # 也可以是数组
  graph::AbstractGraph
  toposort::Vector{Int}
  names::Union{Vector{FT},Nothing} = nothing # 变量名
  index::Vector{CartesianIndex{2}}
  index_rev::Matrix{Int}
  lon::Vector = 1:size(index_rev, 1)
  lat::Vector = 1:size(index_rev, 2)
  nodata::FT = FT(0)
end

function RiverGraph(data::AbstractArray{FT}, g::RiverGraph) where {FT}
  (; graph, toposort, names, index, index_rev, nodata, lon, lat) = g
  RiverGraph{FT}(; graph, toposort, names, data, index, index_rev, nodata, lon, lat)
end


"init a RiverGraph from flowdir matrix"
function RiverGraph(A_fdir::AbstractMatrix{FT}; nodata=FT(0), kw...) where {FT}
  index, index_rev = active_indices(A_fdir, nodata)
  ldd = A_fdir[index]

  graph = graph_flow(ldd, index, pcr_dir)
  toposort = topological_sort_by_dfs(graph)
  RiverGraph(; graph, toposort, data=ldd, nodata, index, index_rev, kw...)
end

function RiverGraph(f::String)
  lon, lat = st_dims(f)
  lat = reverse(lat)
  A = read_flowdir(f)
  RiverGraph(A; lon, lat)
end


"Convert a gridded drainage direction to a directed graph"
function graph_flow(ldd::AbstractVector, inds::AbstractVector, pcr_dir::AbstractVector)
  # prepare a directed graph to be filled
  n = length(inds)
  graph = DiGraph(n)

  # loop over ldd, adding the edge to the downstream node
  for (from_node, from_index) in enumerate(inds)
    ldd_val = ldd[from_node]
    # skip pits to prevent cycles
    ldd_val == 5 && continue
    to_index = from_index + pcr_dir[ldd_val]
    # find the node id of the downstream cell
    to_node = searchsortedfirst(inds, to_index)
    add_edge!(graph, from_node, to_node)
  end
  if is_cyclic(graph)
    error("""One or more cycles detected in flow graph.
        The provided local drainage direction map may be unsound.
        Verify that each active flow cell flows towards a pit.
        """)
  end
  return graph
end



"""
    fillnodata_upstream(g, toposort, data, nodata)

Fill `nodata` upstream cells with the value from the first downstream valid cell, based on
directed acyclic graph `g`, topological order `toposort`, nodata value `nodata` and `data`
containing missing values. Returns filled data `data_out`.
"""
function fillnodata_upstream(g, toposort, data, nodata)
  data_out = copy(data)
  for v in reverse(toposort)  # down- to upstream
    idx_ds = outneighbors(g, v)
    if !isempty(idx_ds)
      if data_out[v] == nodata && data_out[only(idx_ds)] != nodata
        data_out[v] = data_out[only(idx_ds)]
      end
    end
  end
  return data_out
end

fillnodata_upstream(g::RiverGraph, data, nodata) =
  fillnodata_upstream(g.graph, g.toposort, data, nodata)


# 只填充河道的部分
fillnodata_upriver(g::RiverGraph, data, nodata, streamorder; min_sto::Int=4) =
  fillnodata_upriver(g.graph, g.toposort, data, nodata, streamorder; min_sto)

# only for links
function fillnodata_upriver(g::AbstractGraph, toposort,
  links, nodata, streamorder; min_sto::Int=4)

  ids = filter(x -> x > 0, links)
  nodes = findall(x -> x > 0, links)
  find_node(id) = nodes[findfirst(ids .== id)]

  res = copy(links)
  info = []
  # info = Dict{Int,Int}()

  for id_from in reverse(toposort)  # down- to upstream
    streamorder[id_from] < min_sto && continue
    id_to = outneighbors(g, id_from)
    isempty(id_to) && continue

    id_to = only(id_to)
    if res[id_to] != nodata
      if res[id_from] == nodata
        res[id_from] = res[id_to]
      else
        from = res[id_from]
        to = res[id_to]
        # info[from] = to
        push!(info, (; from=find_node(from), to=find_node(to),
          value=from, value_next=to))
      end
    end
  end
  res, DataFrame(info)
end
