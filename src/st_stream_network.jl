import SpatRasters: st_points
# st_points(geoms::Vector{Union{Missing,Shapefile.Point}}) = map(p -> (p.x, p.y), geoms)
st_points(geoms::Vector) = map(p -> (p.x, p.y), geoms)
st_points(shp::DataFrame) = st_points(shp.geometry)


"""
    stream_network(info_node::DataFrame)

# Example
```julia
rg = RiverGraph(f)

min_sto = 4
strord = stream_order(rg)
links = stream_link(rg, strord; min_sto)
river, info_node = fillnodata_upriver(rg, links, strord; min_sto, nodata=0)

stream_network(info_node)
```
"""
function stream_network(info_node::DataFrame)
  vals = info_node |> d -> vcat(d.value, d.value_next)
  n = maximum(vals)
  net_node = DiGraph(n)
  for i = 1:nrow(info_node)
    from, to = info_node[i, [:value, :value_next]]
    add_edge!(net_node, from, to)
  end
  net_node
end


"""
    show_net(net, info_node::DataFrame, info_pour::DataFrame)

 - `info_node`: returned by `fillnodata_upriver` and `flow_path`
 - `info_pour`: with the columns of `site` and `idnex_pit`.
    + `index_pit`: returned by `find_pits`
"""
function show_NetNode(net, info_node::DataFrame, info_pour::DataFrame)
  (; index_pit) = info_pour
  sites = info_pour.site

  for (i, pit) in enumerate(index_pit)
    node = info_node |> d -> d.value_next[d.to.==pit] |> unique
    if !isempty(node)
      r = graph_children(net, only(node))
      println(sites[i], "\t : ", r)
    end
  end
end


function st_stream_network!(rg::RiverGraph, pours; min_sto=5)
  sites = pours.name
  points = st_points(pours)

  # index_pit = find_pits(rg, points) # 这里是往下移动了一个网格
  index_pit = point2index(rg, points)
  info_pour = DataFrame(; site=sites, index_pit)

  strord = stream_order(rg)
  links = stream_link(rg, strord; min_sto)
  # points_link = link2point(rg, links)

  add_links!(links, index_pit)
  ra_basin = fillnodata_upbasin(rg, links; nodata=0)

  river, info_node = fillnodata_upriver(rg, links, strord; min_sto, nodata=0)

  rg.strord .= strord
  rg.links .= links
  rg.river .= river

  net_node = stream_network(info_node)           # 河网结构
  show_NetNode(net_node, info_node, info_pour)
  (; ra_basin, info_node, net_node)
end

export st_stream_network!, stream_network
export show_NetNode
