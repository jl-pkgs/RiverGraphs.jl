
using Graphs, Parameters
import SpatRasters: SpatRaster, st_dims

export find_outlet, graph_children
export read_flowdir

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


# 其他函数调用时，需要知道数据类型，因此此处保留了{FT}
@with_kw mutable struct RiverGraph{FT}
  ngrid::Int
  graph::AbstractGraph

  toposort::Vector{Int} = 1:ngrid
  data::AbstractArray{FT} = zeros(FT, ngrid)# 也可以是数组
  strord::Vector{Int} = zeros(Int, ngrid)
  links::Vector{Int} = zeros(Int, ngrid)
  river::Vector{Int} = zeros(Int, ngrid)
  basins::Vector{Int} = zeros(Int, ngrid)

  names::Union{Vector{FT},Nothing} = nothing # 变量名
  index::Vector{CartesianIndex{2}}
  index_rev::Matrix{Int}
  lon::Vector = 1:size(index_rev, 1)
  lat::Vector = 1:size(index_rev, 2)
  nodata::FT = FT(0)
end

function RiverGraph(data::AbstractArray{FT}, rg::RiverGraph) where {FT}
  (; graph, toposort, names, index, index_rev, nodata, lon, lat) = rg
  ngrid = length(toposort)
  RiverGraph{FT}(; ngrid, graph, toposort, data, names, index, index_rev, nodata, lon, lat)
end


"""
init a RiverGraph from flowdir matrix

- `A`: flowdir
"""
function RiverGraph(A::AbstractMatrix{FT}; nodata=FT(0), kw...) where {FT}
  index, index_rev = active_indices(A, nodata)
  ldd = A[index]

  graph = graph_flow(ldd, index, pcr_dir)
  toposort = topological_sort_by_dfs(graph)
  ngrid = length(toposort)
  RiverGraph(; ngrid, graph, toposort, data=ldd, nodata, index, index_rev, kw...)
end


"""
- `pit`: flowdir is `0`
"""
function read_flowdir(f::String)
  A_gis = read_gdal(f, 1)#[:, end:-1:1] # 修正颠倒的lat
  A = gis2wflow(A_gis)
  nodata = gdal_nodata(f)[1]
  replace!(A, nodata => 5) # nodata as pit
  A
end

# safe: use nodata as pit
function RiverGraph(f::String, points=nothing; safe=true)
  ra = rast(f)
  lon, lat = st_dims(f)
  # A_gis = read_gdal(f, 1)#[:, end:-1:1] # 修正颠倒的lat
  A = gis2wflow(ra.A)
  pit = UInt8(5)

  if !isnothing(points)
    index_pit = point2index(ra, points)
    A[index_pit] .= pit # set pit as 5
  end

  nodata = gdal_nodata(f)[1]
  safe && replace!(A, nodata => pit) # nodata as pit
  RiverGraph(A; lon, lat, nodata)
end

# 重要的教训，流域边界需设置为pit
"Convert a gridded drainage direction to a directed graph"
function graph_flow(ldd::AbstractVector, inds::AbstractVector, pcr_dir::AbstractVector; nodata::Int=0)
  # prepare a directed graph to be filled
  n = length(inds)
  graph = DiGraph(n)

  # loop over ldd, adding the edge to the downstream node
  for (from_node, from_index) in enumerate(inds)
    ldd_val = ldd[from_node]
    # skip pits to prevent cycles
    (ldd_val == 5 || ldd_val == nodata) && continue
    to_index = from_index + pcr_dir[ldd_val]
    # find the node id of the downstream cell
    to_node = searchsortedfirst(inds, to_index)

    # to_node需要在可行的范围内
    if from_node != to_node
      add_edge!(graph, from_node, to_node)
    end
  end

  if is_cyclic(graph)
    error("""One or more cycles detected in flow graph.
        The provided local drainage direction map may be unsound.
        Verify that each active flow cell flows towards a pit.
        """)
  end
  return graph
end

# isvalid_flowdir()

function Base.Matrix(g::RiverGraph, v::AbstractVector{T}, nodata::T=T(0)) where {T}
  reverse_index(v, g.index, g.index_rev; nodata)
end

function Base.Matrix(g::RiverGraph, xs::Tuple, nodata=0)
  map(v -> Matrix(g, v, nodata), xs)
end

Base.Matrix(g::RiverGraph) = Matrix(g, g.data, g.nodata)


function SpatRaster(rg::RiverGraph, x::AbstractVector; nodata=nothing, kw...)
  @assert length(rg.index) == length(x)
  isscalar(nodata) && (nodata = [nodata])

  b = st_bbox(rg.lon, rg.lat)
  A = Matrix(rg, x)
  rast(A, b; nodata, kw...)
end

isscalar(x) = !isa(x, AbstractArray)
isscalar(::Nothing) = false
# export isscalar


## find outlet
_find_outlet(net::SimpleDiGraph) = topological_sort_by_dfs(net)[end] ## 多个节点的时候，该方法会出错

function Base.show(io::IO, net::SimpleDiGraph)
  v = _find_outlet(net)
  println(graph_children(net, v))
end

# 往上追溯一个网格，由于index_pit的流向未设置为0
function find_outlet(net::SimpleDiGraph, toposort, strord; min_sto=2)
  nodes = Vector{Int}()
  orders = Vector{Int}()

  for node in toposort
    node_to = outneighbors(net, node)
    if isempty(node_to) && strord[node] > min_sto
      # node_in = inneighbors(net, node) |> only
      # _node = node_in
      _node = node
      push!(nodes, _node)
      push!(orders, strord[_node])
    end
  end
  nodes, orders
  # DataFrame(; node=nodes, stream_order=)
end

function find_outlet(rg::RiverGraph; min_sto=2)
  (; graph, toposort, strord, links) = rg
  nodes, orders = find_outlet(graph, toposort, strord; min_sto)
  # index = indexin(nodes, toposort)
  DataFrame(; index=nodes, strord=orders, link=links[nodes])
end
