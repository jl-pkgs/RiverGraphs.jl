using SpatRasters: st_cellsize, bbox, st_bbox, bbox_overlap
using SpatRasters: SpatRaster, rast
import SpatRasters: bbox2dims, _zip

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

istrue(x::Bool) = x

find_range(x::BitVector) = findfirst(istrue, x):findlast(istrue, x)

function find_range(mask::BitMatrix)
  inds_x = find_range(sum(mask, dims=2)[:] .!== 0)
  inds_y = find_range(sum(mask, dims=1)[:] .!== 0)
  inds_x, inds_y
end


## 范围取整
function st_range(b::bbox, cellsize)
  cellx, celly = cellsize
  xmin = floor(b.xmin / cellx) * cellx
  xmax = ceil(b.xmax / cellx) * cellx
  ymin = floor(b.ymin / celly) * celly
  ymax = ceil(b.ymax / cellx) * celly
  bbox(; xmin, xmax, ymin, ymax)
end


function st_shrink(ra::SpatRaster{T,N}; nodata=0.0, cellsize_target=nothing) where {T<:Real,N}
  mask = isnan(nodata) ? .!isnan.(ra.A) : ra.A .!== nodata
  lon, lat = st_dims(ra)
  inds, b = st_shrink(mask, lon, lat; cellsize_target)
  ra[inds...]
end


function st_shrink(ra_mask::SpatRaster{Bool}; cellsize_target=nothing)
  lon, lat = st_dims(ra_mask)
  st_shrink(ra_mask.A, lon, lat; cellsize_target)
end

function st_shrink(mask::BitMatrix, lon::AbstractVector, lat::AbstractVector;
  cellsize_target=nothing)
  isnothing(cellsize_target) && (cellsize_target = st_cellsize(lon, lat))
  length(cellsize_target) == 1 && (cellsize_target = (1, 1) .* cellsize_target)

  inds_x, inds_y = find_range(mask)
  reverse_lat = issorted(lat, rev=true)

  box = st_bbox(lon, lat)
  cellsize = st_cellsize(lon, lat)

  _b = st_bbox(lon[inds_x], lat[inds_y])
  b = st_range(_b, cellsize_target)      # 范围稍大一些，要包括所有`b`
  ix, iy = bbox_overlap2(b, box; cellsize, reverse_lat) # ix, iy
  _lon = lon[ix]
  _lat = lat[iy]
  (ix, iy), st_bbox(_lon, _lat)
end
# _lon, _lat, @view data[ix, iy]



function bbox_overlap2(b::bbox, box::bbox; size=nothing, cellsize=nothing, reverse_lat=true, zip=true)
  Lon, Lat = bbox2dims(box; size, cellsize, reverse_lat)

  ilon = findall(b.xmin .< Lon .< b.xmax) |> _zip
  ilat = findall(b.ymin .< Lat .< b.ymax) |> _zip

  if zip
    ilon = _zip(ilon)
    ilat = _zip(ilat)
  end
  # lon, lat = bbox2dims(b; size, cellsize, reverse_lat)
  # @assert length(lon) == length(ilon)
  # @assert length(lat) == length(ilat)
  ilon, ilat
end


export st_shrink
export get_coord, earth_dist, st_length
