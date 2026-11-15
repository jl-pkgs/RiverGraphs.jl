"""
    d8_flowdir(dem; nodata=nothing, cellsize=(1.0, 1.0))

Compute PCRaster/Wflow-style D8 local drainage directions from a DEM.

Direction codes follow `pcr_dir`; `5` denotes a pit/flat and `0` denotes nodata.
The DEM should normally be hydrologically conditioned (breached/filled) before
calling this function. `cellsize` is `(dx, dy)` and only affects slope ranking.
"""
function d8_flowdir(dem::AbstractMatrix; nodata=nothing,
  cellsize::Tuple{<:Real,<:Real}=(1.0, 1.0))
  dx, dy = cellsize
  dx > 0 && dy > 0 || throw(ArgumentError("cellsize must be positive"))

  ldd = fill(UInt8(0), size(dem))
  dist = ntuple(9) do k
    d = pcr_dir[k]
    hypot(abs(d[1]) * dx, abs(d[2]) * dy)
  end

  @inbounds for I in CartesianIndices(dem)
    z = dem[I]
    _invalid_dem_value(z, nodata) && continue

    best_dir = 5
    best_slope = zero(float(z))

    for dir in (1, 2, 3, 4, 6, 7, 8, 9)
      J = I + pcr_dir[dir]
      checkbounds(Bool, dem, J) || continue

      zj = dem[J]
      _invalid_dem_value(zj, nodata) && continue
      drop = z - zj
      drop > 0 || continue

      slope = drop / dist[dir]
      if slope > best_slope
        best_slope = slope
        best_dir = dir
      end
    end

    ldd[I] = UInt8(best_dir)
  end
  ldd
end

_invalid_dem_value(x, ::Nothing) = x isa AbstractFloat && isnan(x)
_invalid_dem_value(x, nodata) = isequal(x, nodata) || (x isa AbstractFloat && isnan(x))


"""
    ldd_code(from, to)

Return the PCRaster/Wflow LDD code for two 8-neighbouring raster cells.
"""
function ldd_code(from::CartesianIndex{2}, to::CartesianIndex{2})
  delta = to - from
  for dir in 1:9
    pcr_dir[dir] == delta && return UInt8(dir)
  end
  throw(ArgumentError("cells $from and $to are not 8-neighbours"))
end


"""
    rasterize_flowpath(x, y, lon, lat)
    rasterize_flowpath(line, lon, lat)

Rasterize an ordered vector river centreline to an ordered, strictly
8-connected grid path. The raster follows RiverGraphs' `(lon, lat)` matrix
layout. Input coordinates and raster axes must use the same CRS.

The line order is preserved. If the source river layer is not already ordered
upstream -> downstream, use `orient_flowpath(...; direction=:dem)` afterwards.
"""
function rasterize_flowpath(x::AbstractVector, y::AbstractVector,
  lon::AbstractVector, lat::AbstractVector)
  length(x) == length(y) || throw(DimensionMismatch("x and y must have equal length"))
  length(x) >= 2 || throw(ArgumentError("a flow path needs at least two vertices"))
  isempty(lon) && throw(ArgumentError("lon axis is empty"))
  isempty(lat) && throw(ArgumentError("lat axis is empty"))

  vertices = Vector{CartesianIndex{2}}(undef, length(x))
  @inbounds for k in eachindex(x, y)
    i = _nearest_axis_index(x[k], lon)
    j = _nearest_axis_index(y[k], lat)
    vertices[k] = CartesianIndex(i, j)
  end

  path = CartesianIndex{2}[]
  for k in 2:length(vertices)
    segment = _grid_line(vertices[k-1], vertices[k])
    for I in segment
      (isempty(path) || path[end] != I) && push!(path, I)
    end
  end

  length(path) >= 2 || throw(ArgumentError("river line collapses to one raster cell"))
  path
end

function rasterize_flowpath(line::AbstractVector{<:Tuple},
  lon::AbstractVector, lat::AbstractVector)
  length(line) >= 2 || throw(ArgumentError("a flow path needs at least two vertices"))
  x = [p[1] for p in line]
  y = [p[2] for p in line]
  rasterize_flowpath(x, y, lon, lat)
end

rasterize_flowpaths(lines::AbstractVector, lon::AbstractVector, lat::AbstractVector) =
  map(line -> rasterize_flowpath(line, lon, lat), lines)


"""
    orient_flowpath(path, dem; direction=:geometry, nodata=nothing)

Orient an 8-connected raster river path.

- `:geometry`: preserve the source line order.
- `:dem` / `:auto`: use the least-squares elevation trend along the full path;
  reverse the path when elevation increases in the source order.

Using the full-path trend is more robust to coarse-DEM noise than comparing only
line endpoints.
"""
function orient_flowpath(path::AbstractVector{CartesianIndex{2}}, dem::AbstractMatrix;
  direction::Symbol=:geometry, nodata=nothing)
  direction === :geometry && return collect(path)
  direction in (:dem, :auto) ||
    throw(ArgumentError("direction must be :geometry, :dem, or :auto"))

  n = 0
  sx = 0.0
  sy = 0.0
  sxx = 0.0
  sxy = 0.0

  @inbounds for (k, I) in enumerate(path)
    checkbounds(Bool, dem, I) || throw(BoundsError(dem, I))
    z = dem[I]
    _invalid_dem_value(z, nodata) && continue
    x = Float64(k)
    y = Float64(z)
    n += 1
    sx += x
    sy += y
    sxx += x * x
    sxy += x * y
  end

  n >= 2 || throw(ArgumentError("at least two valid DEM cells are required to orient a river"))
  denom = n * sxx - sx * sx
  slope = denom == 0 ? 0.0 : (n * sxy - sx * sy) / denom
  slope > 0 ? reverse(collect(path)) : collect(path)
end


"""
    snap_flowpath_junctions(paths; maxdist=1, ambiguous=:error)

Conservatively repair rasterized confluences. Only a path's downstream endpoint
is eligible for snapping, and it may snap only to a cell on another path that
has a downstream continuation. This avoids turning nearby outlet endpoints into
false confluences.

`maxdist` is a Chebyshev distance in raster cells. `0` disables snapping.
Ambiguous equally-near target cells raise an error by default; use
`ambiguous=:skip` to leave such endpoints unchanged.

This is a raster-topology repair heuristic. If the source vector layer carries
an explicit reach-to-reach topology, that topology should take precedence.
"""
function snap_flowpath_junctions(paths::AbstractVector;
  maxdist::Integer=1, ambiguous::Symbol=:error)
  maxdist >= 0 || throw(ArgumentError("maxdist must be non-negative"))
  ambiguous in (:error, :skip) ||
    throw(ArgumentError("ambiguous must be :error or :skip"))

  snapped = map(path -> _compress_flowpath(path), paths)
  maxdist == 0 && return snapped

  path_lengths = length.(snapped)
  owners = Dict{CartesianIndex{2},Vector{Tuple{Int,Int}}}()
  for (j, path) in enumerate(snapped)
    for (k, I) in enumerate(path)
      push!(get!(owners, I, Tuple{Int,Int}[]), (j, k))
    end
  end

  for i in eachindex(snapped)
    path = snapped[i]
    isempty(path) && continue
    endpoint = path[end]

    # Already connected to another reach: nothing to repair.
    shared = false
    for (j, _) in get(owners, endpoint, Tuple{Int,Int}[])
      if j != i
        shared = true
        break
      end
    end
    shared && continue

    best_d2 = typemax(Int)
    targets = Set{CartesianIndex{2}}()

    for di in -maxdist:maxdist, dj in -maxdist:maxdist
      (di == 0 && dj == 0) && continue
      J = endpoint + CartesianIndex(di, dj)
      entries = get(owners, J, Tuple{Int,Int}[])
      isempty(entries) && continue

      eligible = false
      for (j, k) in entries
        # The target reach must continue downstream from the junction cell.
        if j != i && k < path_lengths[j]
          eligible = true
          break
        end
      end
      eligible || continue

      d2 = di * di + dj * dj
      if d2 < best_d2
        best_d2 = d2
        empty!(targets)
        push!(targets, J)
      elseif d2 == best_d2
        push!(targets, J)
      end
    end

    isempty(targets) && continue
    if length(targets) > 1
      ambiguous === :skip && continue
      throw(ArgumentError(
        "ambiguous confluence near $endpoint: equally-near targets $(collect(targets))"))
    end

    target = only(targets)
    bridge = _grid_line(endpoint, target)
    append!(path, @view bridge[2:end])
    snapped[i] = _compress_flowpath(path)
  end

  snapped
end


"""
    condition_river_dem(dem, paths; ...)

Hydrologically condition a DEM around trusted, ordered river paths without
changing terrain away from the river corridor.

The conditioning has three independent parts:

1. `burn_depth` lowers a corridor around river cells, tapered with distance.
2. `bank_drop` guarantees each channel cell is below its adjacent non-channel
   terrain where such terrain exists.
3. `min_slope` enforces a monotonic downstream channel profile over the whole
   river graph. Confluences are handled topologically, so lowering at a tributary
   junction propagates correctly into the downstream reach.

All operations only lower elevations. The returned DEM is `Float64`; the input
is not modified. This is intended for routing, not geomorphic elevation
analysis. Global depression filling/breaching remains a separate preprocessing
step.
"""
function condition_river_dem(dem::AbstractMatrix, paths::AbstractVector;
  nodata=nothing, cellsize::Tuple{<:Real,<:Real}=(1.0, 1.0),
  burn_depth::Real=0.0, burn_width::Integer=0,
  bank_drop::Real=0.01, min_slope::Real=1e-4)

  dx, dy = cellsize
  dx > 0 && dy > 0 || throw(ArgumentError("cellsize must be positive"))
  burn_depth >= 0 || throw(ArgumentError("burn_depth must be non-negative"))
  burn_width >= 0 || throw(ArgumentError("burn_width must be non-negative"))
  bank_drop >= 0 || throw(ArgumentError("bank_drop must be non-negative"))
  min_slope >= 0 || throw(ArgumentError("min_slope must be non-negative"))

  conditioned = Float64.(dem)
  cleanpaths = map(path -> _compress_flowpath(path), paths)
  river = falses(size(dem))

  for path in cleanpaths
    for I in path
      checkbounds(Bool, dem, I) || throw(BoundsError(dem, I))
      _invalid_dem_value(dem[I], nodata) &&
        throw(ArgumentError("river path intersects nodata DEM cell $I"))
      river[I] = true
    end
  end

  # Tapered AGREE-like burn. Each cell is lowered relative to the original DEM,
  # so overlapping river corridors do not accumulate artificial burn depth.
  if burn_depth > 0
    for I in CartesianIndices(river)
      river[I] || continue
      for di in -burn_width:burn_width, dj in -burn_width:burn_width
        J = I + CartesianIndex(di, dj)
        checkbounds(Bool, dem, J) || continue
        _invalid_dem_value(dem[J], nodata) && continue

        r = burn_width == 0 ? 0.0 : hypot(di, dj) / (burn_width + 1)
        r < 1 || continue
        target = Float64(dem[J]) - burn_depth * (1 - r)
        conditioned[J] = min(conditioned[J], target)
      end
    end
  end

  # Make the centreline locally attractive to hillslope D8 without imposing an
  # arbitrary deep trench when the DEM already places the channel correctly.
  if bank_drop > 0
    for I in CartesianIndices(river)
      river[I] || continue
      bank_min = Inf
      for dir in (1, 2, 3, 4, 6, 7, 8, 9)
        J = I + pcr_dir[dir]
        checkbounds(Bool, dem, J) || continue
        river[J] && continue
        _invalid_dem_value(dem[J], nodata) && continue
        bank_min = min(bank_min, conditioned[J])
      end
      isfinite(bank_min) &&
        (conditioned[I] = min(conditioned[I], bank_min - bank_drop))
    end
  end

  min_slope > 0 &&
    _enforce_monotonic_river!(conditioned, cleanpaths; cellsize, min_slope)

  conditioned
end


"""
    force_flowpaths!(ldd, paths; outlet=:keep)

Overwrite D8 directions on ordered river paths while leaving non-river cells
unchanged. Confluences are allowed, but divergent paths are rejected because a
D8 cell can have only one downstream neighbour.

`outlet=:keep` preserves the existing direction at each path endpoint.
`outlet=:pit` sets endpoints that have no downstream continuation in any input
path to PCRaster pit code `5`.
"""
function force_flowpaths!(ldd::AbstractMatrix{<:Integer},
  paths::AbstractVector; outlet::Symbol=:keep)
  outlet in (:keep, :pit) || throw(ArgumentError("outlet must be :keep or :pit"))

  downstream, nodes = _river_downstream(paths)
  endpoints = Set{CartesianIndex{2}}()

  for raw_path in paths
    path = _compress_flowpath(raw_path)
    length(path) >= 2 || continue

    for I in path
      checkbounds(Bool, ldd, I) || throw(BoundsError(ldd, I))
      ldd[I] != 0 || throw(ArgumentError("river path intersects nodata cell $I"))
    end
    push!(endpoints, path[end])
  end

  for (from, to) in downstream
    ldd[from] = ldd_code(from, to)
  end

  if outlet === :pit
    for I in endpoints
      haskey(downstream, I) || (ldd[I] = UInt8(5))
    end
  end

  ldd
end


"""
    river_constrained_flowdir(dem, paths; ...)

Build a D8 flow-direction raster from a coarse DEM while enforcing trusted river
paths.

Information sources are deliberately separated:
- DEM controls hillslope drainage;
- vector hydrography controls channel routing.

By default the river cells are locally conditioned before D8 calculation, then
the final river-cell D8 values are overwritten from the vector paths. Set
`condition=false` to disable DEM conditioning. `junction_radius=1` enables a
conservative one-cell confluence repair after path orientation.

Set `direction=:dem` when path orientation is unknown. `validate=true` builds a
`RiverGraph` once to detect cycles in the final drainage network.
"""
function river_constrained_flowdir(dem::AbstractMatrix,
  paths::AbstractVector{<:AbstractVector{CartesianIndex{2}}};
  nodata=nothing, cellsize::Tuple{<:Real,<:Real}=(1.0, 1.0),
  direction::Symbol=:geometry, junction_radius::Integer=0,
  condition::Bool=true, burn_depth::Real=0.0, burn_width::Integer=0,
  bank_drop::Real=0.01, min_slope::Real=1e-4,
  outlet::Symbol=:keep, validate::Bool=true)

  oriented = map(path -> orient_flowpath(path, dem; direction, nodata), paths)
  connected = junction_radius > 0 ?
    snap_flowpath_junctions(oriented; maxdist=junction_radius) : oriented

  dem_routing = condition ?
    condition_river_dem(dem, connected; nodata, cellsize, burn_depth,
      burn_width, bank_drop, min_slope) : dem

  ldd = d8_flowdir(dem_routing; nodata, cellsize)
  force_flowpaths!(ldd, connected; outlet)

  validate && RiverGraph(ldd; nodata=UInt8(0))
  ldd
end

"""
    river_constrained_flowdir(dem, lines, lon, lat; ...)

Convenience method for vector river lines represented as arrays of `(x, y)`
coordinates. Lines are rasterized to 8-connected paths on the supplied raster
axes before the D8 constraint is applied.
"""
function river_constrained_flowdir(dem::AbstractMatrix, lines::AbstractVector,
  lon::AbstractVector, lat::AbstractVector; kw...)
  size(dem) == (length(lon), length(lat)) ||
    throw(DimensionMismatch("DEM size must equal (length(lon), length(lat))"))
  paths = rasterize_flowpaths(lines, lon, lat)
  river_constrained_flowdir(dem, paths; kw...)
end


function _river_downstream(paths::AbstractVector)
  downstream = Dict{CartesianIndex{2},CartesianIndex{2}}()
  nodes = Set{CartesianIndex{2}}()

  for raw_path in paths
    path = _compress_flowpath(raw_path)
    for I in path
      push!(nodes, I)
    end

    for k in 1:length(path)-1
      from = path[k]
      to = path[k+1]
      ldd_code(from, to) # validates 8-connectivity

      if haskey(downstream, from) && downstream[from] != to
        throw(ArgumentError(
          "divergent river paths at $from: $(downstream[from]) and $to; D8 allows one downstream cell"))
      end
      downstream[from] = to
    end
  end

  downstream, nodes
end

function _enforce_monotonic_river!(dem::AbstractMatrix,
  paths::AbstractVector; cellsize::Tuple{<:Real,<:Real}, min_slope::Real)
  downstream, nodes = _river_downstream(paths)
  isempty(nodes) && return dem

  indegree = Dict{CartesianIndex{2},Int}(I => 0 for I in nodes)
  for (_, to) in downstream
    indegree[to] += 1
  end

  queue = CartesianIndex{2}[]
  sizehint!(queue, length(nodes))
  for I in nodes
    indegree[I] == 0 && push!(queue, I)
  end

  dx, dy = cellsize
  head = 1
  nprocessed = 0
  while head <= length(queue)
    I = queue[head]
    head += 1
    nprocessed += 1

    if haskey(downstream, I)
      J = downstream[I]
      d = J - I
      distance = hypot(abs(d[1]) * dx, abs(d[2]) * dy)
      zmax = dem[I] - min_slope * distance
      dem[J] > zmax && (dem[J] = zmax)

      indegree[J] -= 1
      indegree[J] == 0 && push!(queue, J)
    end
  end

  nprocessed == length(nodes) ||
    error("river paths contain a directed cycle: $nprocessed of $(length(nodes)) cells sorted")
  dem
end

function _compress_flowpath(path)
  out = CartesianIndex{2}[]
  for I in path
    I isa CartesianIndex{2} || throw(ArgumentError("flow paths must contain CartesianIndex{2}"))
    (isempty(out) || out[end] != I) && push!(out, I)
  end
  out
end

function _nearest_axis_index(x, axis::AbstractVector)
  n = length(axis)
  n == 1 && return 1

  asc = axis[end] >= axis[1]
  if asc
    x <= axis[1] && return 1
    x >= axis[end] && return n
  else
    x >= axis[1] && return 1
    x <= axis[end] && return n
  end

  lo = 1
  hi = n
  while hi - lo > 1
    mid = (lo + hi) ÷ 2
    if asc ? axis[mid] <= x : axis[mid] >= x
      lo = mid
    else
      hi = mid
    end
  end
  abs(axis[lo] - x) <= abs(axis[hi] - x) ? lo : hi
end

function _grid_line(a::CartesianIndex{2}, b::CartesianIndex{2})
  x0, y0 = Tuple(a)
  x1, y1 = Tuple(b)
  dx = abs(x1 - x0)
  sx = x0 < x1 ? 1 : -1
  dy = -abs(y1 - y0)
  sy = y0 < y1 ? 1 : -1
  err = dx + dy

  out = CartesianIndex{2}[]
  while true
    push!(out, CartesianIndex(x0, y0))
    x0 == x1 && y0 == y1 && break

    e2 = 2 * err
    if e2 >= dy
      err += dy
      x0 += sx
    end
    if e2 <= dx
      err += dx
      y0 += sy
    end
  end
  out
end


export d8_flowdir, ldd_code
export rasterize_flowpath, rasterize_flowpaths, orient_flowpath
export snap_flowpath_junctions, condition_river_dem
export force_flowpaths!, river_constrained_flowdir
