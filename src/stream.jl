stream_order(g::RiverGraph) = stream_order(g.graph, g.toposort)
stream_link(g::RiverGraph, strord, min_sto) = stream_link(g.graph, g.toposort, strord, min_sto)

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
function stream_link(g, toposort, streamorder, min_sto)
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
river, info_node = fillnodata_upriver(g, links, 0, strord; min_sto)

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
