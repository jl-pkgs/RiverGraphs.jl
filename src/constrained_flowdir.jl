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

Using the full-path trend is more robust to 30 m DEM noise than comparing only
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

  downstream = Dict{CartesianIndex{2},CartesianIndex{2}}()
  endpoints = Set{CartesianIndex{2}}()

  for raw_path in paths
    path = _compress_flowpath(raw_path)
    length(path) >= 2 || continue

    for I in path
      checkbounds(Bool, ldd, I) || throw(BoundsError(ldd, I))
      ldd[I] != 0 || throw(ArgumentError("river path intersects nodata cell $I"))
    end

    for k in 1:length(path)-1
      from = path[k]
      to = path[k+1]
      code = ldd_code(from, to) # also validates 8-connectivity

      if haskey(downstream, from) && downstream[from] != to
        throw(ArgumentError(
          "divergent river paths at $from: $(downstream[from]) and $to; D8 allows one downstream cell"))
      end
      downstream[from] = to
      ldd[from] = code
    end
    push!(endpoints, path[end])
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

Build a D8 flow-direction raster from a DEM, then replace river-cell directions
with high-confidence ordered river paths.

This deliberately separates information sources:
- DEM controls hillslope drainage;
- vector hydrography controls channel routing.

Set `direction=:dem` when path orientation is unknown. `validate=true` builds a
`RiverGraph` once to detect cycles in the final drainage network.
"""
function river_constrained_flowdir(dem::AbstractMatrix,
  paths::AbstractVector{<:AbstractVector{CartesianIndex{2}}};
  nodata=nothing, cellsize::Tuple{<:Real,<:Real}=(1.0, 1.0),
  direction::Symbol=:geometry, outlet::Symbol=:keep, validate::Bool=true)

  ldd = d8_flowdir(dem; nodata, cellsize)
  oriented = map(path -> orient_flowpath(path, dem; direction, nodata), paths)
  force_flowpaths!(ldd, oriented; outlet)

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
export force_flowpaths!, river_constrained_flowdir
