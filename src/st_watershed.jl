import SpatRasters: st_location_exact
export st_watershed

"""
    subbasins(g::RiverGraph; min_sto::Int=4)

# Arguments
- `min_sto`: 过小的sto不认为是stream link

# Return 
- `strord`  : stream order
- `links`   : stream link
- `basinId` : basin Id

# Example
```julia
A = read_flowdir(f) # flowdir, image(A) should looks normal
g = RiverGraph(A)
strord, links, basinId = subbasins(g; min_sto=4)
```
"""
function st_watershed(rg::RiverGraph; level::Int=1, kw...)
  # use river links as pour points
  strord = stream_order(rg)
  links = stream_link(rg, strord; level, kw...)
  fillnodata_upbasin(rg, links; nodata=0)
end


function st_watershed(rg::RiverGraph, points::Vector{Tuple{Float64,Float64}})
  locs = st_location_exact(rg.lon, rg.lat, points)      # 查找index
  index_pit = map(p -> rg.index_rev[p[1], p[2]], locs) # 一维vec对应的顺序

  n_pits = length(index_pit)
  basin = fill(0, length(rg.toposort))
  basin[index_pit] = [1:n_pits;]

  fillnodata_upbasin(rg, basin; nodata=0) # rast
end
