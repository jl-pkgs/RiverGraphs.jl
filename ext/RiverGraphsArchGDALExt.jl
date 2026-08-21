module RiverGraphsArchGDALExt

import ArchGDAL as AG
import RiverGraphs: read_river_lines

"""
    read_river_lines(path::AbstractString; layer=0, strict=true)

Read LineString/MultiLineString features from a GDAL-supported vector file and
return a `Vector{Vector{Tuple{Float64,Float64}}}`. MultiLineString parts are
returned as separate river lines.

Coordinates are not reprojected. The vector data must use the same CRS as the
raster axes passed later to `rasterize_flowpath` / `hydro_enforced_flowdir`.

With `strict=true`, non-line geometries raise an error. With `strict=false`,
unsupported or empty geometries are skipped. GeometryCollection is traversed
recursively, so mixed collections can be handled with `strict=false`.
"""
function read_river_lines(path::AbstractString;
  layer::Integer=0, strict::Bool=true)
  lines = Vector{Vector{Tuple{Float64,Float64}}}()

  AG.read(path) do dataset
    nlayers = AG.nlayer(dataset)
    0 <= layer < nlayers || throw(ArgumentError(
      "layer index $layer is outside 0:$(nlayers - 1) for $path"))

    source = AG.getlayer(dataset, layer)
    for feature in source
      geom = AG.getgeom(feature, 0)
      _append_river_geometry!(lines, geom; strict)
    end
  end

  lines
end

function _append_river_geometry!(lines, geom; strict::Bool)
  name = AG.geomname(geom)
  if ismissing(name)
    strict && throw(ArgumentError("river feature has an empty geometry"))
    return lines
  end

  typename = uppercase(String(name))
  if typename == "LINESTRING" || typename == "LINEARRING"
    n = Int(AG.ngeom(geom))
    if n < 2
      strict && throw(ArgumentError("river LineString has fewer than two points"))
      return lines
    end

    line = Vector{Tuple{Float64,Float64}}(undef, n)
    @inbounds for i in 0:n-1
      point = AG.getpoint(geom, i)
      line[i+1] = (Float64(point[1]), Float64(point[2]))
    end
    push!(lines, line)

  elseif typename == "MULTILINESTRING" || typename == "GEOMETRYCOLLECTION"
    for i in 0:Int(AG.ngeom(geom))-1
      part = AG.getgeom(geom, i)
      _append_river_geometry!(lines, part; strict)
    end

  elseif strict
    throw(ArgumentError(
      "unsupported river geometry $typename; expected LineString or MultiLineString"))
  end

  lines
end

end # module RiverGraphsArchGDALExt
