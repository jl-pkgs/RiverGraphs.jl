import SpatRasters: st_location_exact
export st_watershed

"""
    st_watershed(rg::RiverGraph; level::Int=1, kw...)

# Arguments
- `min_sto`: 过小的sto不认为是stream link

# Return 

# Example
```julia
A = read_flowdir(f) # flowdir, image(A) should looks normal
rg = RiverGraph(A)
ra_basin = st_watershed(rg; min_sto=4)
```
"""
function st_watershed(rg::RiverGraph; level::Int=1, min_sto=nothing, kw...)
  min_sto = get_MinSto(streamorder; level, min_sto)
  # use river links as pour points
  strord = stream_order(rg)
  links = stream_link(rg, strord; min_sto, kw...)
  fillnodata_upbasin(rg, links; nodata=0)
end


function st_watershed(rg::RiverGraph, points::Vector{Tuple{Float64,Float64}})
  index_pit = point2index(rg, points) # 没有往下走一个格点

  n_pits = length(index_pit)
  basin = fill(0, length(rg.toposort))
  basin[index_pit] = [1:n_pits;]
  fillnodata_upbasin(rg, basin; nodata=0) # rast
end
