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
function subbasins(g::RiverGraph; level::Int=1)
  strord = stream_order(g)
  links = stream_link(g, strord; level)
  basinId = fillnodata_upstream(g, links, 0)
  strord, links, basinId
end

export subbasins
