"""
    priority_flood_dem(dem, seeds; ...)
    priority_flood_dem(dem, paths; ...)

Create a depression-resolved routing DEM with the Priority-Flood algorithm.
Trusted river cells can be supplied as `seeds`, or indirectly as ordered river
`paths`. Seed elevations are never changed; all other valid cells are visited
outward from the seeds and raised only when necessary to preserve a drainage
path back to an already processed cell.

`min_slope > 0` adds a tiny positive gradient across filled flats so a later D8
calculation has an unambiguous lower neighbour. `boundary_outlets=true` also
uses valid cells on the outer raster edge as outlets. Set it to `false` for a
basin/domain that should drain only to the supplied river network.

This function is intended to create a routing surface. It does not represent a
geomorphically corrected DEM.
"""
function priority_flood_dem(dem::AbstractMatrix,
  seeds::AbstractVector{CartesianIndex{2}};
  nodata=nothing, cellsize::Tuple{<:Real,<:Real}=(1.0, 1.0),
  min_slope::Real=1e-4, boundary_outlets::Bool=true)

  dx, dy = cellsize
  dx > 0 && dy > 0 || throw(ArgumentError("cellsize must be positive"))
  min_slope >= 0 || throw(ArgumentError("min_slope must be non-negative"))

  T = eltype(dem) <: AbstractFloat ? eltype(dem) : Float64
  flooded = T.(dem)
  visited = falses(size(dem))

  heap_z = Float64[]
  heap_i = Int[]
  LI = LinearIndices(dem)
  CI = CartesianIndices(dem)

  function add_seed!(I::CartesianIndex{2})
    checkbounds(Bool, dem, I) || throw(BoundsError(dem, I))
    _invalid_dem_value(dem[I], nodata) && return false
    visited[I] && return false

    visited[I] = true
    _minheap_push!(heap_z, heap_i, Float64(flooded[I]), LI[I])
    true
  end

  for I in seeds
    add_seed!(I)
  end

  if boundary_outlets
    nx, ny = size(dem)
    for i in 1:nx
      add_seed!(CartesianIndex(i, 1))
      ny > 1 && add_seed!(CartesianIndex(i, ny))
    end
    for j in 1:ny
      add_seed!(CartesianIndex(1, j))
      nx > 1 && add_seed!(CartesianIndex(nx, j))
    end
  end

  isempty(heap_z) &&
    throw(ArgumentError("priority flood requires at least one valid seed or boundary outlet"))

  while !isempty(heap_z)
    z_from, linear_from = _minheap_pop!(heap_z, heap_i)
    I = CI[linear_from]

    for dir in (1, 2, 3, 4, 6, 7, 8, 9)
      J = I + pcr_dir[dir]
      checkbounds(Bool, dem, J) || continue
      visited[J] && continue
      _invalid_dem_value(dem[J], nodata) && continue

      visited[J] = true
      d = pcr_dir[dir]
      distance = hypot(abs(d[1]) * dx, abs(d[2]) * dy)
      z_required = z_from + min_slope * distance
      z_new = max(Float64(flooded[J]), z_required)

      # Preserve a representable positive gradient for Float32 routing DEMs.
      value = T(z_new)
      if min_slope > 0 && Float64(value) <= z_from
        value = nextfloat(T(z_from))
      end
      flooded[J] = value
      _minheap_push!(heap_z, heap_i, Float64(value), LI[J])
    end
  end

  # With river-only outlets, disconnected valid components are a data/topology
  # problem and should not silently become pits.
  for I in CartesianIndices(dem)
    _invalid_dem_value(dem[I], nodata) && continue
    visited[I] || throw(ArgumentError(
      "valid DEM component containing $I is disconnected from all priority-flood outlets"))
  end

  flooded
end

function priority_flood_dem(dem::AbstractMatrix,
  paths::AbstractVector{<:AbstractVector{CartesianIndex{2}}}; kw...)
  _, nodes = _river_downstream(paths)
  seeds = collect(nodes)
  priority_flood_dem(dem, seeds; kw...)
end


# Minimal binary min-heap. Keeping keys and linear raster indices in separate
# primitive arrays avoids a DataStructures dependency and is substantially more
# memory-efficient than heap entries containing CartesianIndex tuples on very
# large rasters.
function _minheap_push!(keys::Vector{Float64}, values::Vector{Int},
  key::Float64, value::Int)
  push!(keys, key)
  push!(values, value)
  i = length(keys)

  while i > 1
    parent = i >>> 1
    keys[parent] <= keys[i] && break
    keys[parent], keys[i] = keys[i], keys[parent]
    values[parent], values[i] = values[i], values[parent]
    i = parent
  end
  nothing
end

function _minheap_pop!(keys::Vector{Float64}, values::Vector{Int})
  isempty(keys) && throw(ArgumentError("cannot pop an empty heap"))

  key = keys[1]
  value = values[1]
  last_key = pop!(keys)
  last_value = pop!(values)

  if !isempty(keys)
    keys[1] = last_key
    values[1] = last_value
    i = 1
    n = length(keys)

    while true
      left = i << 1
      left > n && break
      right = left + 1
      child = right <= n && keys[right] < keys[left] ? right : left
      keys[i] <= keys[child] && break

      keys[i], keys[child] = keys[child], keys[i]
      values[i], values[child] = values[child], values[i]
      i = child
    end
  end

  key, value
end


export priority_flood_dem
