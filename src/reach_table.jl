function _node_values(rg::RiverGraph, data, name)
  values = data isa SpatRaster ? data.A : data
  if ndims(values) == 2 && size(values) == size(rg.index_rev)
    return values[rg.index]
  end
  length(values) == rg.ngrid ||
    throw(DimensionMismatch("$name must match the flow-direction grid"))
  vec(values)
end

function _distance_m(rg::RiverGraph, from::Int, to::Int, lonlat::Bool)
  i1, j1 = Tuple(rg.index[from])
  i2, j2 = Tuple(rg.index[to])
  x1, y1 = rg.lon[i1], rg.lat[j1]
  x2, y2 = rg.lon[i2], rg.lat[j2]
  lonlat ? 1000 * earth_dist(x1, y1, x2, y2) : hypot(x2 - x1, y2 - y1)
end

function _cell_area_m2(rg::RiverGraph, node::Int, lonlat::Bool)
  i, j = Tuple(rg.index[node])
  dx, dy = _cellspan(rg.lon, i), _cellspan(rg.lat, j)
  lonlat || return dx * dy

  radius = 6_371_008.8
  φ1 = deg2rad(rg.lat[j] - dy / 2)
  φ2 = deg2rad(rg.lat[j] + dy / 2)
  radius^2 * deg2rad(dx) * abs(sin(φ2) - sin(φ1))
end

function _cell_areas(rg::RiverGraph, cell_area_m2, lonlat::Bool)
  area = if isnothing(cell_area_m2)
    [_cell_area_m2(rg, node, lonlat) for node in 1:rg.ngrid]
  elseif cell_area_m2 isa Number
    fill(Float64(cell_area_m2), rg.ngrid)
  else
    Float64.(_node_values(rg, cell_area_m2, "cell_area_m2"))
  end
  all(x -> isfinite(x) && x > 0, area) ||
    throw(ArgumentError("cell_area_m2 must contain positive finite values"))
  area
end

"""
    delineate_reaches(rg; river_mask=nothing, min_sto=2,
        elevation=nothing, cell_area_m2=nothing)

按源头、每个汇流点和出口划分河段，生成 Reach/HRU Table、河段栅格和
HRU 编号栅格。每个 HRU 是汇入对应 Reach 的局部汇水区。未提供
`river_mask` 时，以 `stream_order >= min_sto` 提取河道；Shapefile 应先
栅格化至流向网格。`elevation` 单位为 m。
经纬度坐标按球面距离计算，投影坐标默认单位为 m。
"""
function delineate_reaches(rg::RiverGraph; river_mask=nothing, min_sto::Int=2,
  elevation=nothing, cell_area_m2=nothing)
  min_sto > 0 || throw(ArgumentError("min_sto must be positive"))
  graph, order = rg.graph, rg.toposort

  river = if isnothing(river_mask)
    stream_order(rg) .>= min_sto
  else
    values = _node_values(rg, river_mask, "river_mask")
    nodata = river_mask isa SpatRaster ? first(river_mask.nodata) : nothing
    map(x -> !ismissing(x) && !isequal(x, nodata) && !iszero(x), values)
  end
  any(river) || throw(ArgumentError("river network is empty"))

  upstream_count = zeros(Int, rg.ngrid)
  for node in order
    river[node] || continue
    upstream_count[node] = count(upstream -> river[upstream], inneighbors(graph, node))
  end

  starts = [node for node in order if river[node] && upstream_count[node] != 1]
  paths = Vector{Vector{Int}}()
  start_reach = zeros(Int, rg.ngrid)
  for start in starts
    start_reach[start] = length(paths) + 1
    path = Int[start]
    node = start
    while true
      downstream = outneighbors(graph, node)
      isempty(downstream) && break
      node = only(downstream)
      river[node] || break
      push!(path, node)
      upstream_count[node] != 1 && break
    end
    push!(paths, path)
  end

  nreach = length(paths)
  river_reach = zeros(Int, rg.ngrid)
  downstream_id = zeros(Int, nreach)
  for (id, path) in pairs(paths)
    endpoint = last(path)
    next_id = endpoint == first(path) ? 0 : start_reach[endpoint]
    downstream_id[id] = next_id
    stop = next_id == 0 ? length(path) : length(path) - 1
    river_reach[path[1:stop]] .= id
  end
  all(>(0), river_reach[river]) || error("river mask is not connected by the flow direction")

  lonlat = all(x -> -180 <= x <= 180, rg.lon) && all(y -> -90 <= y <= 90, rg.lat)
  lengths = [sum(_distance_m(rg, path[i - 1], path[i], lonlat)
    for i in 2:length(path); init=0.0) for path in paths]
  slopes = fill(NaN, nreach)
  if !isnothing(elevation)
    values = _node_values(rg, elevation, "elevation")
    nodata = elevation isa SpatRaster ? first(elevation.nodata) : nothing
    any(node -> river[node] &&
      (ismissing(values[node]) || isequal(values[node], nodata)), eachindex(river)) &&
      throw(ArgumentError("elevation contains nodata on the river network"))
    z = Float64.(values)
    for (id, path) in pairs(paths)
      lengths[id] > 0 && (slopes[id] = max(0.0, (z[first(path)] - z[last(path)]) / lengths[id]))
    end
  end

  hru_id = fillnodata_upbasin(graph, order, river_reach; nodata=0)
  cell_area = _cell_areas(rg, cell_area_m2, lonlat)
  areas = zeros(nreach)
  for node in eachindex(hru_id)
    id = hru_id[node]
    id == 0 || (areas[id] += cell_area[node])
  end

  reaches = DataFrame(
    segId=1:nreach,
    downSegId=downstream_id,
    length=lengths,
    slope=slopes,
    from_node=first.(paths),
    to_node=last.(paths),
    n_node=length.(paths),
    index=paths,
  )
  hrus = DataFrame(HRUid=1:nreach, hruSegId=1:nreach, area=areas)
  return (;
    reaches,
    hrus,
    river_reach=SpatRaster(rg, river_reach; nodata=0),
    hru_id=SpatRaster(rg, hru_id; nodata=0),
  )
end

function _grid_runoff(rg::RiverGraph, runoff_depth)
  shape = size(rg.index_rev)
  if ndims(runoff_depth) == 3 && size(runoff_depth)[1:2] == shape
    runoff = Matrix{Float64}(undef, rg.ngrid, size(runoff_depth, 3))
    for t in axes(runoff_depth, 3), node in 1:rg.ngrid
      i, j = Tuple(rg.index[node])
      runoff[node, t] = runoff_depth[i, j, t]
    end
  elseif ndims(runoff_depth) == 2 && size(runoff_depth) == shape
    runoff = reshape(Float64.(runoff_depth[rg.index]), rg.ngrid, 1)
  elseif ndims(runoff_depth) == 2 && size(runoff_depth, 1) == rg.ngrid
    runoff = Float64.(runoff_depth)
  elseif ndims(runoff_depth) == 1 && length(runoff_depth) == rg.ngrid
    runoff = reshape(Float64.(runoff_depth), rg.ngrid, 1)
  else
    throw(DimensionMismatch("runoff_depth must be [lon, lat, time] or [ngrid, time]"))
  end
  all(isfinite, runoff) || throw(ArgumentError("runoff_depth contains non-finite values"))
  all(>=(0), runoff) || throw(ArgumentError("runoff_depth contains negative values"))
  runoff
end

"""
    aggregate_hru_runoff(rg, runoff_depth; kwargs...)

将网格径流深率按面积加权到局部 HRU，返回 mizuRoute 所需的 Reach Table、
HRU Table 和 `HRU × time` 径流矩阵。网格维度为 `[lon, lat, time]`，
也可传入按 `rg.index` 排列的 `[ngrid, time]`；径流单位保持不变。
"""
function aggregate_hru_runoff(rg::RiverGraph, runoff_depth; river_mask=nothing,
  min_sto::Int=2, elevation=nothing, cell_area_m2=nothing)
  result = delineate_reaches(rg; river_mask, min_sto, elevation, cell_area_m2)
  runoff = _grid_runoff(rg, runoff_depth)
  hru_id = result.hru_id.A[rg.index]
  all(>(0), hru_id) ||
    throw(ArgumentError("some active cells do not drain to the extracted river network"))
  lonlat = all(x -> -180 <= x <= 180, rg.lon) && all(y -> -90 <= y <= 90, rg.lat)
  cell_area = _cell_areas(rg, cell_area_m2, lonlat)

  hru_runoff = zeros(nrow(result.hrus), size(runoff, 2))
  for node in 1:rg.ngrid
    id = hru_id[node]
    id == 0 && continue
    for t in axes(runoff, 2)
      hru_runoff[id, t] += runoff[node, t] * cell_area[node]
    end
  end
  hru_runoff ./= result.hrus.area

  return (;
    reaches=result.reaches,
    hrus=result.hrus,
    runoff=hru_runoff,
    river_reach=result.river_reach,
    hru_id=result.hru_id,
  )
end

export delineate_reaches, aggregate_hru_runoff
