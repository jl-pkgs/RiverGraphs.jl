export getInfo_links

stream_order(g::RiverGraph) = stream_order(g.graph, g.toposort)
stream_link(g::RiverGraph, strord; level) = stream_link(g.graph, g.toposort, strord; level)

"""
    stream_order(g, toposort)

Return the Strahler streamorder based on directed acyclic graph `g` and topological order
`toposort`.
"""
function stream_order(g, toposort)
  n = length(toposort)
  strord = fill(1, n)
  for v in toposort
    inds_up = inneighbors(g, v)
    if !isempty(inds_up)
      sto_up = strord[inds_up]
      if length(findall(x -> x == maximum(sto_up), sto_up)) > 1
        strord[v] = (maximum(sto_up) + 1)
      else
        strord[v] = maximum(sto_up)
      end
    end
  end
  return strord
end

"""
    stream_link(g, streamorder, toposort, min_sto)

Return stream_link with a unique id starting at 1, based on a minimum streamorder `min_sto`,
directed acyclic graph `g` and topological order `toposort`.
"""
function stream_link(g, toposort, streamorder; level::Int=2)
  min_sto = max(maximum(streamorder) - level, 1)

  n = length(toposort)
  links = fill(0, n)

  i = 1
  for v in toposort
    streamorder[v] < min_sto && continue

    ds_nodes = outneighbors(g, v)
    if !isempty(ds_nodes)
      if streamorder[v] != streamorder[only(ds_nodes)]
        links[v] = i
        i += 1
      end
    else
      # also set pits (without a downstream node)
      links[v] = i
      i += 1
    end
  end
  return links
end


function getInfo_links(g::RiverGraph, links_2d)
  con = links_2d .!= 0
  xs, ys = get_coord(con)
  vals = filter(x -> x != 0, links_2d)
  inds = g.index_rev[con][:]
  DataFrame(; x=xs, y=ys,
    lon=g.lon[xs], lat=g.lat[ys],
    link=vals, index=inds)
end

"""
    stream_network(info_node::DataFrame)

# Example
```julia
A = read_flowdir(path_flowdir_GuanShan)
g = RiverGraph(A)

min_sto = 4
@time strord, links, basinId = subbasins(g; min_sto)
strord_2d, links_2d, basinId_2d =
  Matrix(g, strord, -1), Matrix(g, links), Matrix(g, basinId)
river, info_node = fillnodata_upriver(g, links, strord; min_sto, nodata=0)

flow_path(g, info_node, strord; min_sto)
```
"""
function stream_network(info_node::DataFrame)
  vals = info_node |> d -> vcat(d.value, d.value_next)
  n = maximum(vals)
  graph = DiGraph(n)
  for i = 1:nrow(info_node)
    from, to = info_node[i, [:value, :value_next]]
    add_edge!(graph, from, to)
  end
  graph
end
