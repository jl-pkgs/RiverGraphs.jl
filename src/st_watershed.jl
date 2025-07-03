import SpatRasters: st_location_exact
export st_watershed, st_subbasins, write_subbasins

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


function st_subbasins(ra_basin::SpatRaster, info_node::DataFrame, ids::AbstractVector;
  cellsize_target=0.1)

  A2 = deepcopy(ra_basin.A)
  for I in eachindex(A2)
    if !(A2[I] in ids)
      A2[I] = 0
    end
  end
  r = SpatRaster(A2, ra_basin)
  r2 = st_shrink(r; cellsize_target)

  inds = findall(i -> i in ids, info_node.value)
  info_node = info_node[inds, :]
  (; basin=r2, info_node)
end


"""
    write_subbasins(rg::RiverGraph, info_node::DataFrame, pours::DataFrame)

- `pours`: with columns of `name` and `geometry`
"""
function write_subbasins(rg::RiverGraph, info_node::DataFrame, pours::DataFrame;
  outdir="./OUTPUT")

  sites = pours.name
  points = st_points(pours)

  # index = find_pits(rg, points)
  index = point2index(rg, points)
  ids = rg.links[index]
  @show ids
  # info_pour = DataFrame(; link=ids, site=sites, geometry=points)
  net_node = stream_network(info_node)
  ra_basin = SpatRaster(rg, rg.basins)

  N = length(ids)
  for i in 1:N
    node = ids[i]
    _r = graph_children(net_node, node)
    _ids = unlist(_r)
    isscalar(_ids) && (_ids = [_ids])

    _basin, _info_node = st_subbasins(ra_basin, info_node, _ids)
    # _net = stream_network(_info_node)

    site = sites[i]
    println("[$i] $site\t: $_r")

    f_tif = @sprintf("%s/%02d_%s_basinId.tif", outdir, i, site)
    f_csv = @sprintf("%s/%02d_%s_info_node.csv", outdir, i, site)
    write_gdal(_basin, f_tif; nodata=0)
    fwrite(_info_node, f_csv)
  end
end
