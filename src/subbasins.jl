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
strord_2d, links_2d, basinId_2d = subbasins(g; min_sto=4)
```
"""
function subbasins(g::RiverGraph; min_sto::Int=4)
  strord = stream_order(g)
  links = stream_link(g, strord, min_sto)
  basinId = fillnodata_upstream(g, links, 0)

  strord, links, basinId
  # Matrix(g, (strord, links, basinId))
  # Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
end

# function subbasins(A_fdir::AbstractMatrix{FT};
#   min_sto::Int=4, nodata::FT=FT(0)) where {FT<:Real}

#   inds, inds_rev = active_indices(A_fdir, nodata)
#   g = graph_flow(A_fdir[inds], inds, pcr_dir)
#   toposort = topological_sort_by_dfs(g)

#   strord = stream_order(g, toposort)
#   strord_2d = reverse_index(strord, inds, inds_rev; nodata=-1)

#   links = stream_link(g, toposort, strord, min_sto)
#   links_2d = reverse_index(links, inds, inds_rev; nodata=0)

#   basinId = fillnodata_upstream(g, toposort, links, 0)
#   basinId_2d = reverse_index(basinId, inds, inds_rev; nodata=0)

#   strord_2d, links_2d, basinId_2d
# end


export subbasins
