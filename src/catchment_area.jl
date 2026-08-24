export upstream_nodes, catchment_area, catchment_stats

"""返回出水口及其全部上游节点。"""
function upstream_nodes(rg::RiverGraph, outlet::Integer)
  checkbounds(rg.toposort, outlet)
  nodes = Int[outlet]
  seen = falses(rg.ngrid)
  seen[outlet] = true

  for i in eachindex(nodes)
    for node in inneighbors(rg.graph, nodes[i])
      seen[node] && continue
      seen[node] = true
      push!(nodes, node)
    end
  end
  nodes
end

_cellspan(v, i) = length(v) == 1 ? error("坐标轴至少需要两个点") :
  i == firstindex(v) ? abs(v[i + 1] - v[i]) :
  i == lastindex(v) ? abs(v[i] - v[i - 1]) : abs(v[i + 1] - v[i - 1]) / 2

function _cell_area_km2(rg::RiverGraph, node::Int)
  i, j = Tuple(rg.index[node])
  dx, dy = _cellspan(rg.lon, i), _cellspan(rg.lat, j)
  lonlat = all(x -> -180 <= x <= 180, rg.lon) && all(y -> -90 <= y <= 90, rg.lat)
  lonlat || return dx * dy / 1e6

  radius = 6_371_008.8
  φ1 = deg2rad(rg.lat[j] - dy / 2)
  φ2 = deg2rad(rg.lat[j] + dy / 2)
  radius^2 * deg2rad(dx) * abs(sin(φ2) - sin(φ1)) / 1e6
end

"""
    catchment_area(rg, outlet; cell_area_km2=nothing)

计算出水口上游集水面积（km²）。`cell_area_km2` 可覆盖坐标推算值。
经纬度坐标按球面网格计算；投影坐标默认单位为米。
"""
function catchment_area(rg::RiverGraph, outlet::Integer; cell_area_km2=nothing)
  nodes = upstream_nodes(rg, outlet)
  isnothing(cell_area_km2) ? sum(node -> _cell_area_km2(rg, node), nodes) :
    length(nodes) * cell_area_km2
end

function catchment_area(rg::RiverGraph, point::Tuple{<:Real,<:Real}; kw...)
  outlet = only(point2index(rg, [Float64.(point)]))
  outlet > 0 || error("出水口不在有效流域内：$point")
  catchment_area(rg, outlet; kw...)
end

"""返回出水口节点、上游节点数和集水面积。"""
function catchment_stats(rg::RiverGraph, outlet::Integer; kw...)
  nodes = upstream_nodes(rg, outlet)
  (; outlet, upstream_nodes=length(nodes), area_km2=catchment_area(rg, outlet; kw...))
end

function catchment_stats(rg::RiverGraph, point::Tuple{<:Real,<:Real}; kw...)
  outlet = only(point2index(rg, [Float64.(point)]))
  outlet > 0 || error("出水口不在有效流域内：$point")
  catchment_stats(rg, outlet; kw...)
end
