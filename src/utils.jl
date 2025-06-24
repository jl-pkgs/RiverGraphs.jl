"Map from PCRaster LDD value to a CartesianIndex"
const pcr_dir = [
  CartesianIndex(-1, -1),  # 1
  CartesianIndex(0, -1),  # 2
  CartesianIndex(1, -1),  # 3
  CartesianIndex(-1, 0),  # 4
  CartesianIndex(0, 0),  # 5
  CartesianIndex(1, 0),  # 6
  CartesianIndex(-1, 1),  # 7
  CartesianIndex(0, 1),  # 8
  CartesianIndex(1, 1),  # 9
]


"""
    active_indices(subcatch_2d, nodata)

Takes a 2D array of the subcatchments. And derive forward and reverse indices.

1: Get a list of `CartesianIndex{2}`` that are active, based on a nodata value.
These map from the 1D internal domain to the 2D external domain.

2: Make a reverse index, a `Matrix{Int}``, which maps from the 2D external domain to
the 1D internal domain, providing an Int which can be used as a linear index. Values of 0
represent inactive cells.
"""
function active_indices(A::AbstractMatrix, nodata)
  all_inds = CartesianIndices(size(A))
  indices = filter(i -> !isequal(A[i], nodata), all_inds) # 不为NA的全部收纳

  reverse_indices = zeros(Int, size(A))
  for (i, I) in enumerate(indices)
    reverse_indices[I] = i
  end
  return indices, reverse_indices
end


function reverse_index(v::AbstractVector{T}, index, index_rev; nodata::T=T(0)) where {T}
  R = zeros(T, size(index_rev)) .+ nodata
  @inbounds for (i, I) in enumerate(index)
    R[I] = v[i]
  end
  R
end

function Base.Matrix(g::RiverGraph, v::AbstractVector{T}, nodata::T=T(0)) where {T}
  reverse_index(v, g.index, g.index_rev; nodata)
end

function Base.Matrix(g::RiverGraph, xs::Tuple, nodata=0)
  map(v -> Matrix(g, v, nodata), xs)
end

Base.Matrix(g::RiverGraph) = Matrix(g, g.data, g.nodata)


function get_coord(inds::Vector{CartesianIndex{2}})
  xs = map(p -> p[1], inds)
  ys = map(p -> p[2], inds)
  xs, ys
end
get_coord(lgl::BitArray) = get_coord(findall(lgl))

# check manually
function meshgrid(x, y)
  X = repeat(x, 1, length(y))
  Y = repeat(y', length(x), 1)
  X, Y
end

function earth_dist(lon1::T, lat1::T, lon2::T, lat2::T; R=6371.0) where {T<:Real}
  φ1 = deg2rad(lat1)
  φ2 = deg2rad(lat2)
  Δφ = deg2rad(lat2 - lat1)
  Δλ = deg2rad(lon2 - lon1)

  a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
  c = 2 * atan(sqrt(a), sqrt(1 - a))
  return R * c
end

function earth_dist(p1::Tuple{T,T}, p2::Tuple{T,T}; R=6371.0) where {T<:Real}
  lat1, lon1 = p1
  lat2, lon2 = p2
  earth_dist(lat1, lon1, lat2, lon2; R)
end

function st_length(lon::AbstractVector, lat::AbstractVector)
  dist = 0.0
  for i = 2:lastindex(lon)
    p1 = lon[i-1], lat[i-1]
    p2 = lon[i], lat[i]
    dist += earth_dist(p1, p2)
  end
  dist # in Km 
end

# function st_length(LON::AbstractMatrix, LAT::AbstractMatrix, inds)
#   st_length(LON[inds], LAT[inds])
# end


export get_coord, reverse_index, earth_dist, st_length
